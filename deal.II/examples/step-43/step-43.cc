/* $Id$ */
/* Author: Chih-Che Chueh, University of Victoria, 2010 */
/*         Wolfgang Bangerth, Texas A&M University, 2010 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2010 by Chih-Che Chueh and the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/utilities.h>
#include <base/function.h>
#include <base/tensor_function.h>

#include <lac/full_matrix.h>
#include <lac/solver_gmres.h>
#include <lac/solver_cg.h>
#include <lac/block_sparsity_pattern.h>
#include <lac/constraint_matrix.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_tools.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>

#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>

#include <numerics/vectors.h>
#include <numerics/data_out.h>
#include <numerics/solution_transfer.h>

#include <lac/trilinos_sparse_matrix.h>
#include <lac/trilinos_block_sparse_matrix.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_block_vector.h>
#include <lac/trilinos_precondition.h>

#include <fstream>
#include <sstream>

using namespace dealii;

namespace LinearSolvers
{
  template <class Matrix, class Preconditioner>
  class InverseMatrix : public Subscriptor
  {
    public:
      InverseMatrix (const Matrix         &m,
                     const Preconditioner &preconditioner);


      template <typename VectorType>
      void vmult (VectorType       &dst,
                  const VectorType &src) const;

    private:
      const SmartPointer<const Matrix> matrix;
      const Preconditioner &preconditioner;
  };


  template <class Matrix, class Preconditioner>
  InverseMatrix<Matrix,Preconditioner>::
  InverseMatrix (const Matrix &m,
                 const Preconditioner &preconditioner)
                 :
                 matrix (&m),
                 preconditioner (preconditioner)
  {}



  template <class Matrix, class Preconditioner>
  template <typename VectorType>
  void
  InverseMatrix<Matrix,Preconditioner>::
  vmult (VectorType       &dst,
         const VectorType &src) const
  {
    SolverControl solver_control (src.size(), 1e-7*src.l2_norm());
    SolverCG<VectorType> cg (solver_control);

    dst = 0;

    try
      {
       cg.solve (*matrix, dst, src, preconditioner);
      }
    catch (std::exception &e)
      {
       Assert (false, ExcMessage(e.what()));
      }
  }

  template <class PreconditionerA, class PreconditionerMp>
  class BlockSchurPreconditioner : public Subscriptor
  {
    public:
      BlockSchurPreconditioner (
        const TrilinosWrappers::BlockSparseMatrix     &S,
        const InverseMatrix<TrilinosWrappers::SparseMatrix,
                            PreconditionerMp>         &Mpinv,
        const PreconditionerA                         &Apreconditioner);

      void vmult (TrilinosWrappers::BlockVector       &dst,
                  const TrilinosWrappers::BlockVector &src) const;

    private:
      const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> darcy_matrix;
      const SmartPointer<const InverseMatrix<TrilinosWrappers::SparseMatrix,
                                             PreconditionerMp > > m_inverse;
      const PreconditionerA &a_preconditioner;

      mutable TrilinosWrappers::Vector tmp;
  };



  template <class PreconditionerA, class PreconditionerMp>
  BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::
  BlockSchurPreconditioner(const TrilinosWrappers::BlockSparseMatrix  &S,
                           const InverseMatrix<TrilinosWrappers::SparseMatrix,
                                               PreconditionerMp>      &Mpinv,
                           const PreconditionerA                      &Apreconditioner)
                  :
                  darcy_matrix            (&S),
                  m_inverse               (&Mpinv),
                  a_preconditioner        (Apreconditioner),
                  tmp                     (darcy_matrix->block(1,1).m())
  {}


  template <class PreconditionerA, class PreconditionerMp>
  void BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::vmult (
     TrilinosWrappers::BlockVector       &dst,
     const TrilinosWrappers::BlockVector &src) const
  {
    a_preconditioner.vmult (dst.block(0), src.block(0));
    darcy_matrix->block(1,0).residual(tmp, dst.block(0), src.block(1));
    tmp *= -1;
    m_inverse->vmult (dst.block(1), tmp);
  }
}


template <int dim>
class TwoPhaseFlowProblem
{
  public:
    TwoPhaseFlowProblem (const unsigned int degree);
    void run ();

  private:
    void setup_dofs ();
    void assemble_darcy_preconditioner ();
    void build_darcy_preconditioner ();
    void assemble_darcy_system ();
    void assemble_saturation_system ();
    void assemble_saturation_matrix ();
    void assemble_saturation_rhs ();
    void assemble_saturation_rhs_cell_term (const FEValues<dim>             &saturation_fe_values,
                                            const FEValues<dim>             &darcy_fe_values,
                                            const std::vector<unsigned int> &local_dof_indices,
                                            const double                     global_u_infty,
                                            const double                     global_S_variation,
                                            const double                     global_Omega_diameter);
    void assemble_saturation_rhs_boundary_term (const FEFaceValues<dim>             &saturation_fe_face_values,
                                                const FEFaceValues<dim>             &darcy_fe_face_values,
                                                const std::vector<unsigned int>     &local_dof_indices);
    double get_maximal_velocity () const;
    std::pair<double,double> get_extrapolated_saturation_range () const;
    void solve ();
    bool determine_whether_to_solve_pressure_velocity_part () const;
    void compute_refinement_indicators (Vector<double> &indicator) const;
    void refine_grid (const Vector<double> &indicator);
    void project_back_saturation ();
    void output_results () const;

    static
    double
    compute_viscosity(const std::vector<double>          &old_saturation,
                      const std::vector<double>          &old_old_saturation,
                      const std::vector<Tensor<1,dim> >  &old_saturation_grads,
                      const std::vector<Tensor<1,dim> >  &old_old_saturation_grads,
                      const std::vector<Vector<double> > &present_darcy_values,
                      const double                        global_u_infty,
                      const double                        global_S_variation,
                      const double                        global_Omega_diameter,
                      const double                        cell_diameter,
                      const double                        old_time_step,
                      const double                        viscosity);


    const unsigned int degree;

    Triangulation<dim>                   triangulation;

    const unsigned int                   darcy_degree;
    FESystem<dim>                        darcy_fe;
    DoFHandler<dim>                      darcy_dof_handler;
    ConstraintMatrix                     darcy_constraints;

    ConstraintMatrix                     darcy_preconditioner_constraints;

    TrilinosWrappers::BlockSparseMatrix  darcy_matrix;
    TrilinosWrappers::BlockSparseMatrix  darcy_preconditioner_matrix;

    TrilinosWrappers::BlockVector        darcy_solution;
    TrilinosWrappers::BlockVector        darcy_rhs;

    TrilinosWrappers::BlockVector        nth_darcy_solution_after_solving_pressure_part;
    TrilinosWrappers::BlockVector        n_minus_oneth_darcy_solution_after_solving_pressure_part;

    const unsigned int                   saturation_degree;
    FE_Q<dim>                            saturation_fe;
    DoFHandler<dim>                      saturation_dof_handler;
    ConstraintMatrix                     saturation_constraints;

    TrilinosWrappers::SparseMatrix       saturation_matrix;

    TrilinosWrappers::Vector             predictor_saturation_solution;
    TrilinosWrappers::Vector             saturation_solution;
    TrilinosWrappers::Vector             old_saturation_solution;
    TrilinosWrappers::Vector             old_old_saturation_solution;
    TrilinosWrappers::Vector             saturation_rhs;

    TrilinosWrappers::Vector             nth_saturation_solution_after_solving_pressure_part;

    const unsigned int        n_refinement_steps;
    bool                      solve_pressure_velocity_part;
    bool                      previous_solve_pressure_velocity_part;

    const double              saturation_level;
    const double              saturation_value;

    double n_minus_oneth_time_step;
    double cumulative_nth_time_step;

    double time_step;
    double old_time_step;
    unsigned int timestep_number;
    double viscosity;

    std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionIC> Amg_preconditioner;
    std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionIC>  Mp_preconditioner;

    bool rebuild_saturation_matrix;
};


template <int dim>
class PressureRightHandSide : public Function<dim>
{
  public:
    PressureRightHandSide () : Function<dim>(1) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
};



template <int dim>
double
PressureRightHandSide<dim>::value (const Point<dim>  &/*p*/,
                                   const unsigned int /*component*/) const
{
  return 0;
}


template <int dim>
class PressureBoundaryValues : public Function<dim>
{
  public:
    PressureBoundaryValues () : Function<dim>(1) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
};


template <int dim>
double
PressureBoundaryValues<dim>::value (const Point<dim>  &p,
                                    const unsigned int /*component*/) const
{
  return 1-p[0];
}


template <int dim>
class SaturationBoundaryValues : public Function<dim>
{
  public:
    SaturationBoundaryValues () : Function<dim>(1) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
};



template <int dim>
double
SaturationBoundaryValues<dim>::value (const Point<dim> &p,
                                      const unsigned int /*component*/) const
{
  if (p[0] == 0)
    return 1;
  else
    return 0;
}


template <int dim>
class SaturationInitialValues : public Function<dim>
{
  public:
    SaturationInitialValues () : Function<dim>(1) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;

};


template <int dim>
double
SaturationInitialValues<dim>::value (const Point<dim>  &/*p*/,
                                     const unsigned int /*component*/) const
{
  return 0;
}


template <int dim>
void
SaturationInitialValues<dim>::vector_value (const Point<dim> &p,
                                            Vector<double>   &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = SaturationInitialValues<dim>::value (p,c);
}


namespace SingleCurvingCrack
{
  template <int dim>
  class KInverse : public TensorFunction<2,dim>
  {
    public:
      KInverse ()
		      :
		      TensorFunction<2,dim> ()
	{}

      virtual void value_list (const std::vector<Point<dim> > &points,
                               std::vector<Tensor<2,dim> >    &values) const;
  };


  template <int dim>
  void
  KInverse<dim>::value_list (const std::vector<Point<dim> > &points,
                             std::vector<Tensor<2,dim> >    &values) const
  {
    Assert (points.size() == values.size(),
            ExcDimensionMismatch (points.size(), values.size()));

    for (unsigned int p=0; p<points.size(); ++p)
      {
        values[p].clear ();

        const double distance_to_flowline
          = std::fabs(points[p][1]-0.5-0.1*std::sin(10*points[p][0]));

        const double permeability = std::max(std::exp(-(distance_to_flowline*
                                                        distance_to_flowline)
                                                      / (0.1 * 0.1)),
                                             0.01);

        for (unsigned int d=0; d<dim; ++d)
          values[p][d][d] = 1./permeability;
      }
  }
}


namespace RandomMedium
{
  template <int dim>
  class KInverse : public TensorFunction<2,dim>
  {
    public:
      KInverse ()
		      :
		      TensorFunction<2,dim> ()
	{}

      virtual void value_list (const std::vector<Point<dim> > &points,
                               std::vector<Tensor<2,dim> >    &values) const;

    private:
      static std::vector<Point<dim> > centers;

      static std::vector<Point<dim> > get_centers ();
  };



  template <int dim>
  std::vector<Point<dim> >
  KInverse<dim>::centers = KInverse<dim>::get_centers();


  template <int dim>
  std::vector<Point<dim> >
  KInverse<dim>::get_centers ()
  {
    const unsigned int N = (dim == 2 ?
                            40 :
                            (dim == 3 ?
                             100 :
                             throw ExcNotImplemented()));

    std::vector<Point<dim> > centers_list (N);
    for (unsigned int i=0; i<N; ++i)
      for (unsigned int d=0; d<dim; ++d)
        centers_list[i][d] = static_cast<double>(rand())/RAND_MAX;

    return centers_list;
  }



  template <int dim>
  void
  KInverse<dim>::value_list (const std::vector<Point<dim> > &points,
                             std::vector<Tensor<2,dim> >    &values) const
  {
    Assert (points.size() == values.size(),
            ExcDimensionMismatch (points.size(), values.size()));

    for (unsigned int p=0; p<points.size(); ++p)
      {
        values[p].clear ();

        double permeability = 0;
        for (unsigned int i=0; i<centers.size(); ++i)
          permeability += std::exp(-(points[p]-centers[i]).square()
                                   / (0.05 * 0.05));

        const double normalized_permeability
          = std::min (std::max(permeability, 0.01), 4.);

        for (unsigned int d=0; d<dim; ++d)
          values[p][d][d] = 1./normalized_permeability;
      }
  }
}


double mobility_inverse (const double S,
                         const double viscosity)
{
  return 1.0 /(1.0/viscosity * S * S + (1-S) * (1-S));
}

double f_saturation (const double S,
                     const double viscosity)
{
  return S*S /( S * S +viscosity * (1-S) * (1-S));
}

double get_fractional_flow_derivative (const double S,
                                       const double viscosity)
{
  const double temp = ( S * S + viscosity * (1-S) * (1-S) );

  const double numerator   =  2.0 * S * temp
                              -
                              S * S *
                              ( 2.0 * S - 2.0 * viscosity * (1-S) );

  const double denomerator =  std::pow(temp, 2.0 );

  return numerator / denomerator;
}

template <int dim>
TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem (const unsigned int degree)
                :
                degree (degree),
                darcy_degree (degree),
                darcy_fe (FE_Q<dim>(darcy_degree+1), dim,
                          FE_Q<dim>(darcy_degree), 1),
                darcy_dof_handler (triangulation),

                saturation_degree (degree),
                saturation_fe (saturation_degree),
                saturation_dof_handler (triangulation),

                n_refinement_steps (4),
                solve_pressure_velocity_part (false),
                previous_solve_pressure_velocity_part (false),

                saturation_level (2),
                saturation_value (0.5),

                time_step (0),
                old_time_step (0),
                viscosity (0.2),

                rebuild_saturation_matrix (true)
{}


template <int dim>
void TwoPhaseFlowProblem<dim>::setup_dofs ()
{
  std::vector<unsigned int> darcy_block_component (dim+1,0);
  darcy_block_component[dim] = 1;
  {
    darcy_dof_handler.distribute_dofs (darcy_fe);
    DoFRenumbering::Cuthill_McKee (darcy_dof_handler);
    DoFRenumbering::component_wise (darcy_dof_handler, darcy_block_component);

    darcy_constraints.clear ();
    DoFTools::make_hanging_node_constraints (darcy_dof_handler, darcy_constraints);
    darcy_constraints.close ();
  }
  {
    saturation_dof_handler.distribute_dofs (saturation_fe);

    saturation_constraints.clear ();
    DoFTools::make_hanging_node_constraints (saturation_dof_handler, saturation_constraints);
    saturation_constraints.close ();
  }
  {
    darcy_preconditioner_constraints.clear ();

    std::vector<bool> component_mask (dim+1, false);
    component_mask[dim] = true;


    DoFTools::make_hanging_node_constraints (darcy_dof_handler, darcy_preconditioner_constraints);
    DoFTools::make_zero_boundary_constraints (darcy_dof_handler, darcy_preconditioner_constraints, component_mask);

    darcy_preconditioner_constraints.close ();
  }


  std::vector<unsigned int> darcy_dofs_per_block (2);
  DoFTools::count_dofs_per_block (darcy_dof_handler, darcy_dofs_per_block, darcy_block_component);
  const unsigned int n_u = darcy_dofs_per_block[0],
                     n_p = darcy_dofs_per_block[1],
                     n_s = saturation_dof_handler.n_dofs();

  std::cout << "Number of active cells: "
            << triangulation.n_active_cells()
            << " (on "
            << triangulation.n_levels()
            << " levels)"
            << std::endl
            << "Number of degrees of freedom: "
            << n_u + n_p + n_s
            << " (" << n_u << '+' << n_p << '+'<< n_s <<')'
            << std::endl
            << std::endl;

  {
    darcy_matrix.clear ();

    BlockCompressedSimpleSparsityPattern csp (2,2);

    csp.block(0,0).reinit (n_u, n_u);
    csp.block(0,1).reinit (n_u, n_p);
    csp.block(1,0).reinit (n_p, n_u);
    csp.block(1,1).reinit (n_p, n_p);

    csp.collect_sizes ();

    Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);

    for (unsigned int c=0; c<dim+1; ++c)
      for (unsigned int d=0; d<dim+1; ++d)
        if (! ((c==dim) && (d==dim)))
          coupling[c][d] = DoFTools::always;
        else
          coupling[c][d] = DoFTools::none;


    DoFTools::make_sparsity_pattern (darcy_dof_handler, coupling, csp,
                                     darcy_constraints, false);

    darcy_matrix.reinit (csp);
  }

  {
    Amg_preconditioner.reset ();
    Mp_preconditioner.reset ();
    darcy_preconditioner_matrix.clear ();

    BlockCompressedSimpleSparsityPattern csp (2,2);

    csp.block(0,0).reinit (n_u, n_u);
    csp.block(0,1).reinit (n_u, n_p);
    csp.block(1,0).reinit (n_p, n_u);
    csp.block(1,1).reinit (n_p, n_p);

    csp.collect_sizes ();

    Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
    for (unsigned int c=0; c<dim+1; ++c)
      for (unsigned int d=0; d<dim+1; ++d)
        if (c == d)
          coupling[c][d] = DoFTools::always;
        else
          coupling[c][d] = DoFTools::none;

    DoFTools::make_sparsity_pattern (darcy_dof_handler, coupling, csp,
                                     darcy_constraints, false);

    darcy_preconditioner_matrix.reinit (csp);
  }


  {
    saturation_matrix.clear ();

    CompressedSimpleSparsityPattern csp (n_s, n_s);

    DoFTools::make_sparsity_pattern (saturation_dof_handler, csp,
                                     saturation_constraints, false);


    saturation_matrix.reinit (csp);
  }

  darcy_solution.reinit (2);
  darcy_solution.block(0).reinit (n_u);
  darcy_solution.block(1).reinit (n_p);
  darcy_solution.collect_sizes ();

  nth_darcy_solution_after_solving_pressure_part.reinit (2);
  nth_darcy_solution_after_solving_pressure_part.block(0).reinit (n_u);
  nth_darcy_solution_after_solving_pressure_part.block(1).reinit (n_p);
  nth_darcy_solution_after_solving_pressure_part.collect_sizes ();

  n_minus_oneth_darcy_solution_after_solving_pressure_part.reinit (2);
  n_minus_oneth_darcy_solution_after_solving_pressure_part.block(0).reinit (n_u);
  n_minus_oneth_darcy_solution_after_solving_pressure_part.block(1).reinit (n_p);
  n_minus_oneth_darcy_solution_after_solving_pressure_part.collect_sizes ();

  darcy_rhs.reinit (2);
  darcy_rhs.block(0).reinit (n_u);
  darcy_rhs.block(1).reinit (n_p);
  darcy_rhs.collect_sizes ();

  predictor_saturation_solution.reinit (n_s);
  saturation_solution.reinit (n_s);
  old_saturation_solution.reinit (n_s);
  old_old_saturation_solution.reinit (n_s);

  nth_saturation_solution_after_solving_pressure_part.reinit (n_s);

  saturation_rhs.reinit (n_s);
}

template <int dim>
void
TwoPhaseFlowProblem<dim>::assemble_darcy_preconditioner ()
{
  std::cout << "   Rebuilding darcy preconditioner..." << std::endl;

  darcy_preconditioner_matrix = 0;

  const QGauss<dim> quadrature_formula(darcy_degree+2);
  FEValues<dim>     darcy_fe_values (darcy_fe, quadrature_formula,
                                     update_JxW_values |
                                     update_values |
                                     update_gradients |
                                     update_quadrature_points);
  FEValues<dim> saturation_fe_values (saturation_fe, quadrature_formula, 
                                      update_values);
 
  const unsigned int   dofs_per_cell   = darcy_fe.dofs_per_cell;
  const unsigned int   n_q_points      = quadrature_formula.size();

  const RandomMedium::KInverse<dim> k_inverse;   
//  const SingleCurvingCrack::KInverse<dim> k_inverse;

  std::vector<Tensor<2,dim> >       k_inverse_values (n_q_points);
  Tensor<2,dim>                     k_value;

  std::vector<double>               old_saturation_values (n_q_points);
 
  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  std::vector<Tensor<1,dim> > phi_u   (dofs_per_cell); 
  std::vector<Tensor<1,dim> > grad_phi_p (dofs_per_cell);
 
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);
 
  typename DoFHandler<dim>::active_cell_iterator
    cell = darcy_dof_handler.begin_active(),
    endc = darcy_dof_handler.end();
  typename DoFHandler<dim>::active_cell_iterator
    saturation_cell = saturation_dof_handler.begin_active();

  for (; cell!=endc; ++cell, ++saturation_cell)
    {
      darcy_fe_values.reinit (cell);
      saturation_fe_values.reinit (saturation_cell);

      local_matrix = 0;

      saturation_fe_values.get_function_values (old_saturation_solution, old_saturation_values);

      k_inverse.value_list (darcy_fe_values.get_quadrature_points(),
                            k_inverse_values);
 
      for (unsigned int q=0; q<n_q_points; ++q)
        {
          const double old_s = old_saturation_values[q];
          const double mobility = 1.0 / mobility_inverse(old_s,viscosity);

          k_value.clear ();
          for (unsigned int d=0; d<dim; d++)
            k_value[d][d] = 1.0 / k_inverse_values[q][d][d];           

          for (unsigned int k=0; k<dofs_per_cell; ++k)
            {
              phi_u[k]       = darcy_fe_values[velocities].value (k,q);
              grad_phi_p[k]  = darcy_fe_values[pressure].gradient (k,q);
            }
 
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                local_matrix(i,j) += (k_inverse_values[q] * mobility_inverse(old_s,viscosity) *
                                      phi_u[i] * phi_u[j]
                                      +
                                      k_value * mobility *
                                      grad_phi_p[i] * grad_phi_p[j])
                                     * darcy_fe_values.JxW(q);
              }
        }
 
      cell->get_dof_indices (local_dof_indices);
      darcy_preconditioner_constraints.distribute_local_to_global (local_matrix,
                                                                   local_dof_indices,
                                                                   darcy_preconditioner_matrix);
    }
}

template <int dim>
void
TwoPhaseFlowProblem<dim>::build_darcy_preconditioner ()
{
  assemble_darcy_preconditioner ();

  Amg_preconditioner = std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionIC>
                       (new TrilinosWrappers::PreconditionIC());
  Amg_preconditioner->initialize(darcy_preconditioner_matrix.block(0,0));

  Mp_preconditioner = std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionIC>
                      (new TrilinosWrappers::PreconditionIC());
  Mp_preconditioner->initialize(darcy_preconditioner_matrix.block(1,1));

}

template <int dim>
void TwoPhaseFlowProblem<dim>::assemble_darcy_system ()
{
  darcy_matrix = 0;
  darcy_rhs    = 0;

  QGauss<dim>   quadrature_formula(darcy_degree+2);
  QGauss<dim-1> face_quadrature_formula(darcy_degree+2);

  FEValues<dim> darcy_fe_values (darcy_fe, quadrature_formula,
                                 update_values    | update_gradients |
                                 update_quadrature_points  | update_JxW_values);

  FEValues<dim> saturation_fe_values (saturation_fe, quadrature_formula,
                                      update_values);

  FEFaceValues<dim> darcy_fe_face_values (darcy_fe, face_quadrature_formula,
                                          update_values    | update_normal_vectors |
                                          update_quadrature_points  | update_JxW_values);

  const unsigned int   dofs_per_cell   = darcy_fe.dofs_per_cell;

  const unsigned int   n_q_points      = quadrature_formula.size();
  const unsigned int   n_face_q_points = face_quadrature_formula.size();

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  const PressureRightHandSide<dim>  pressure_right_hand_side;
  const PressureBoundaryValues<dim> pressure_boundary_values;
  const RandomMedium::KInverse<dim> k_inverse;
//  const SingleCurvingCrack::KInverse<dim> k_inverse;

  std::vector<double>               pressure_rhs_values (n_q_points);
  std::vector<double>               boundary_values (n_face_q_points);
  std::vector<Tensor<2,dim> >       k_inverse_values (n_q_points);

  std::vector<double>               old_saturation_values (n_q_points);

  std::vector<Tensor<1,dim> > phi_u (dofs_per_cell);
  std::vector<double>         div_phi_u (dofs_per_cell);
  std::vector<double>         phi_p (dofs_per_cell);

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  typename DoFHandler<dim>::active_cell_iterator
    cell = darcy_dof_handler.begin_active(),
    endc = darcy_dof_handler.end();
  typename DoFHandler<dim>::active_cell_iterator
    saturation_cell = saturation_dof_handler.begin_active();

  for (; cell!=endc; ++cell, ++saturation_cell)
    {
      darcy_fe_values.reinit (cell);
      saturation_fe_values.reinit (saturation_cell);

      local_matrix = 0;
      local_rhs = 0;

      saturation_fe_values.get_function_values (old_saturation_solution, old_saturation_values);

      pressure_right_hand_side.value_list (darcy_fe_values.get_quadrature_points(),
                                           pressure_rhs_values);
      k_inverse.value_list (darcy_fe_values.get_quadrature_points(),
                            k_inverse_values);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int k=0; k<dofs_per_cell; ++k)
            {
              phi_u[k]     = darcy_fe_values[velocities].value (k,q);
              div_phi_u[k] = darcy_fe_values[velocities].divergence (k,q);
              phi_p[k]     = darcy_fe_values[pressure].value (k,q);
            }
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              const double old_s = old_saturation_values[q];
              for (unsigned int j=0; j<=i; ++j)
                {
                  local_matrix(i,j) += (phi_u[i] * k_inverse_values[q] *
                                        mobility_inverse(old_s,viscosity) * phi_u[j]
                                        - div_phi_u[i] * phi_p[j]
                                        - phi_p[i] * div_phi_u[j])
                                       * darcy_fe_values.JxW(q);
                }

              local_rhs(i) += (-phi_p[i] * pressure_rhs_values[q])*
                               darcy_fe_values.JxW(q);
            }
         }

      for (unsigned int face_no=0;
           face_no<GeometryInfo<dim>::faces_per_cell;
           ++face_no)
        if (cell->at_boundary(face_no))
          {
            darcy_fe_face_values.reinit (cell, face_no);

            pressure_boundary_values
              .value_list (darcy_fe_face_values.get_quadrature_points(),
                           boundary_values);

            for (unsigned int q=0; q<n_face_q_points; ++q)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  const Tensor<1,dim>
                    phi_i_u = darcy_fe_face_values[velocities].value (i, q);

                  local_rhs(i) += -(phi_i_u *
                                    darcy_fe_face_values.normal_vector(q) *
                                    boundary_values[q] *
                                    darcy_fe_face_values.JxW(q));
                }
          }

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=i+1; j<dofs_per_cell; ++j)
          local_matrix(i,j) = local_matrix(j,i);

      cell->get_dof_indices (local_dof_indices);

      darcy_constraints.distribute_local_to_global (local_matrix,
                                                    local_rhs,
                                                    local_dof_indices,
                                                    darcy_matrix,
                                                    darcy_rhs);

    }
}

template <int dim>
void TwoPhaseFlowProblem<dim>::assemble_saturation_system ()
{
  if ( rebuild_saturation_matrix == true )
    {
      saturation_matrix = 0;
      assemble_saturation_matrix ();
    }

  saturation_rhs = 0;
  assemble_saturation_rhs ();
}

template <int dim>
void TwoPhaseFlowProblem<dim>::assemble_saturation_matrix ()
{
  QGauss<dim> quadrature_formula(saturation_degree+2);

  FEValues<dim> saturation_fe_values (saturation_fe, quadrature_formula,
                                      update_values | update_JxW_values);

  const unsigned int dofs_per_cell = saturation_fe.dofs_per_cell;

  const unsigned int n_q_points = quadrature_formula.size();

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
    cell = saturation_dof_handler.begin_active(),
    endc = saturation_dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      saturation_fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs    = 0;

      for (unsigned int q=0; q<n_q_points; ++q)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const double phi_i_s = saturation_fe_values.shape_value (i,q);
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const double phi_j_s = saturation_fe_values.shape_value (j,q);
                local_matrix(i,j) += phi_i_s * phi_j_s * saturation_fe_values.JxW(q);
              }
          }
      cell->get_dof_indices (local_dof_indices);

      saturation_constraints.distribute_local_to_global (local_matrix,
                                                         local_dof_indices,
                                                         saturation_matrix);

    }
}


template <int dim>
void TwoPhaseFlowProblem<dim>::assemble_saturation_rhs ()
{
  QGauss<dim>   quadrature_formula(saturation_degree+2);
  QGauss<dim-1> face_quadrature_formula(saturation_degree+2);

  FEValues<dim> saturation_fe_values                          (saturation_fe, quadrature_formula,
                                                               update_values    | update_gradients |
                                                               update_quadrature_points  | update_JxW_values);
  FEValues<dim> darcy_fe_values                               (darcy_fe, quadrature_formula,
                                                               update_values);
  FEFaceValues<dim> saturation_fe_face_values                 (saturation_fe, face_quadrature_formula,
                                                               update_values    | update_normal_vectors |
                                                               update_quadrature_points  | update_JxW_values);
  FEFaceValues<dim> darcy_fe_face_values                      (darcy_fe, face_quadrature_formula,
                                                               update_values);
  FEFaceValues<dim> saturation_fe_face_values_neighbor        (saturation_fe, face_quadrature_formula,
                                                               update_values);

  const unsigned int dofs_per_cell = saturation_dof_handler.get_fe().dofs_per_cell;
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  const double global_u_infty = get_maximal_velocity ();
  const std::pair<double,double>
    global_S_range = get_extrapolated_saturation_range ();
  const double global_S_variasion = global_S_range.second - global_S_range.first;
  const double global_Omega_diameter = GridTools::diameter (triangulation);

  typename DoFHandler<dim>::active_cell_iterator
    cell = saturation_dof_handler.begin_active(),
    endc = saturation_dof_handler.end();
  typename DoFHandler<dim>::active_cell_iterator
    darcy_cell = darcy_dof_handler.begin_active();
  for (; cell!=endc; ++cell, ++darcy_cell)
    {
      saturation_fe_values.reinit (cell);
      darcy_fe_values.reinit (darcy_cell);

      cell->get_dof_indices (local_dof_indices);

      assemble_saturation_rhs_cell_term(saturation_fe_values,
                                        darcy_fe_values,
                                        local_dof_indices,
                                        global_u_infty,
                                        global_S_variasion,
                                        global_Omega_diameter);

      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
           ++face_no)
        {

          if (cell->at_boundary(face_no))
            {
              darcy_fe_face_values.reinit (darcy_cell, face_no);
              saturation_fe_face_values.reinit (cell, face_no);
              assemble_saturation_rhs_boundary_term (saturation_fe_face_values,
                                                     darcy_fe_face_values,
                                                     local_dof_indices);
            }
        }
    }
}

template <int dim>
void
TwoPhaseFlowProblem<dim>::
assemble_saturation_rhs_cell_term (const FEValues<dim>             &saturation_fe_values,
                                   const FEValues<dim>             &darcy_fe_values,
                                   const std::vector<unsigned int> &local_dof_indices,
                                   const double                     global_u_infty,
                                   const double                     global_S_variation,
                                   const double                     global_Omega_diameter)
{
  const unsigned int dofs_per_cell = saturation_fe_values.dofs_per_cell;
  const unsigned int n_q_points    = saturation_fe_values.n_quadrature_points;

  Vector<double> local_rhs (dofs_per_cell);

  std::vector<double>          old_saturation_solution_values(n_q_points);
  std::vector<double>          old_old_saturation_solution_values(n_q_points);
  std::vector<Tensor<1,dim> >  old_grad_saturation_solution_values(n_q_points);
  std::vector<Tensor<1,dim> >  old_old_grad_saturation_solution_values(n_q_points);
  std::vector<Vector<double> > present_darcy_solution_values(n_q_points, Vector<double>(dim+1));

  saturation_fe_values.get_function_values (old_saturation_solution, old_saturation_solution_values);
  saturation_fe_values.get_function_values (old_old_saturation_solution, old_old_saturation_solution_values);
  saturation_fe_values.get_function_grads (old_saturation_solution, old_grad_saturation_solution_values);
  saturation_fe_values.get_function_grads (old_old_saturation_solution, old_old_grad_saturation_solution_values);
  darcy_fe_values.get_function_values (darcy_solution, present_darcy_solution_values);

  const double nu
    = compute_viscosity (old_saturation_solution_values,
                         old_old_saturation_solution_values,
                         old_grad_saturation_solution_values,
                         old_old_grad_saturation_solution_values,
                         present_darcy_solution_values,
                         global_u_infty,
                         global_S_variation,
                         global_Omega_diameter,
                         saturation_fe_values.get_cell()->diameter(),
                         old_time_step,
                         viscosity);

  for (unsigned int q=0; q<n_q_points; ++q)
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
        const double old_s = old_saturation_solution_values[q];
        Tensor<1,dim> present_u;
        for (unsigned int d=0; d<dim; ++d)
          present_u[d] = present_darcy_solution_values[q](d);

        const double        phi_i_s      = saturation_fe_values.shape_value (i, q);
        const Tensor<1,dim> grad_phi_i_s = saturation_fe_values.shape_grad (i, q);

        local_rhs(i) += (time_step *
                         f_saturation(old_s,viscosity) *
                         present_u *
                         grad_phi_i_s
                         -
                         time_step *
                         nu *
                         old_grad_saturation_solution_values[q] * grad_phi_i_s
                         +
                         old_s * phi_i_s)
                         *
                         saturation_fe_values.JxW(q);
      }

  saturation_constraints.distribute_local_to_global (local_rhs,
                                                     local_dof_indices,
                                                     saturation_rhs);
}

template <int dim>
void
TwoPhaseFlowProblem<dim>::
assemble_saturation_rhs_boundary_term (const FEFaceValues<dim>             &saturation_fe_face_values,
                                       const FEFaceValues<dim>             &darcy_fe_face_values,
                                       const std::vector<unsigned int>     &local_dof_indices)
{
  const unsigned int dofs_per_cell      = saturation_fe_face_values.dofs_per_cell;
  const unsigned int n_face_q_points    = saturation_fe_face_values.n_quadrature_points;

  Vector<double> local_rhs (dofs_per_cell);

  std::vector<double>          old_saturation_solution_values_face(n_face_q_points);
  std::vector<Vector<double> > present_darcy_solution_values_face(n_face_q_points, Vector<double>(dim+1));
  std::vector<double>          neighbor_saturation (n_face_q_points);

  saturation_fe_face_values.get_function_values (old_saturation_solution, old_saturation_solution_values_face);
  darcy_fe_face_values.get_function_values (darcy_solution, present_darcy_solution_values_face);

  SaturationBoundaryValues<dim> saturation_boundary_values;
  saturation_boundary_values
    .value_list (saturation_fe_face_values.get_quadrature_points(),
                 neighbor_saturation);

  for (unsigned int q=0; q<n_face_q_points; ++q)
    {
      Tensor<1,dim> present_u_face;
      for (unsigned int d=0; d<dim; ++d)
        present_u_face[d] = present_darcy_solution_values_face[q](d);

      const double normal_flux = present_u_face *
                                 saturation_fe_face_values.normal_vector(q);

      const bool is_outflow_q_point = (normal_flux >= 0);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        local_rhs(i) -= time_step *
                        normal_flux *
                        f_saturation((is_outflow_q_point == true
                        ?
                        old_saturation_solution_values_face[q]
                        :
                        neighbor_saturation[q]),
                        viscosity) *
                        saturation_fe_face_values.shape_value (i,q) *
                        saturation_fe_face_values.JxW(q);
    }
  saturation_constraints.distribute_local_to_global (local_rhs,
                                                     local_dof_indices,
                                                     saturation_rhs);
}


template <int dim>
void TwoPhaseFlowProblem<dim>::solve ()
{
  solve_pressure_velocity_part = determine_whether_to_solve_pressure_velocity_part ();

  if ( timestep_number <= 3 || solve_pressure_velocity_part == true )
    {
      std::cout << "   Solving darcy system (pressure-velocity part)..." << std::endl;

      assemble_darcy_system ();
      build_darcy_preconditioner ();

      {
        const LinearSolvers::InverseMatrix<TrilinosWrappers::SparseMatrix,
                                           TrilinosWrappers::PreconditionIC>
          mp_inverse (darcy_preconditioner_matrix.block(1,1), *Mp_preconditioner);

        const LinearSolvers::BlockSchurPreconditioner<TrilinosWrappers::PreconditionIC,
                                                      TrilinosWrappers::PreconditionIC>
          preconditioner (darcy_matrix, mp_inverse, *Amg_preconditioner);

        SolverControl solver_control (darcy_matrix.m(),
                                      1e-6*darcy_rhs.l2_norm());

        SolverGMRES<TrilinosWrappers::BlockVector>
          gmres (solver_control,
                 SolverGMRES<TrilinosWrappers::BlockVector >::AdditionalData(100));

        for (unsigned int i=0; i<darcy_solution.size(); ++i)
          if (darcy_constraints.is_constrained(i))
            darcy_solution(i) = 0;

        gmres.solve(darcy_matrix, darcy_solution, darcy_rhs, preconditioner);

        darcy_constraints.distribute (darcy_solution);

        std::cout << "     "
                  << solver_control.last_step()
                  << " GMRES iterations for darcy system (pressure-velocity part)."
                  << std::endl;

      }

      {
        n_minus_oneth_darcy_solution_after_solving_pressure_part = nth_darcy_solution_after_solving_pressure_part;
        nth_darcy_solution_after_solving_pressure_part = darcy_solution;

        nth_saturation_solution_after_solving_pressure_part = saturation_solution;
      }
    }
  else
    {
      darcy_solution.block(0) = nth_darcy_solution_after_solving_pressure_part.block(0);
      darcy_solution.block(0).sadd (2.0, -1.0, n_minus_oneth_darcy_solution_after_solving_pressure_part.block(0) );

      double extrapolated_time_step = GridTools::minimal_cell_diameter(triangulation) /
                                      get_maximal_velocity() / 8.0;

      double local_cumulative_time_step = cumulative_nth_time_step + extrapolated_time_step;
      double coef_1 = local_cumulative_time_step / n_minus_oneth_time_step;
      double coef_2 = ( 1.0 + coef_1 );

      TrilinosWrappers::Vector tmp (darcy_solution.block(0).size());
      tmp = nth_darcy_solution_after_solving_pressure_part.block(0);

      tmp.sadd (coef_2, -coef_1, n_minus_oneth_darcy_solution_after_solving_pressure_part.block(0) );

      darcy_solution.block(0).sadd (0.5, 0.5, tmp);
    }


  old_time_step = time_step;
  time_step = GridTools::minimal_cell_diameter(triangulation) /
              get_maximal_velocity() / 8.0;

  if ( timestep_number <= 3 || ( solve_pressure_velocity_part == true && previous_solve_pressure_velocity_part == true ) )
    {
      n_minus_oneth_time_step = time_step;
      cumulative_nth_time_step = 0.0;
    }
  else if ( solve_pressure_velocity_part == true && previous_solve_pressure_velocity_part == false )
    {
      n_minus_oneth_time_step = cumulative_nth_time_step;
      cumulative_nth_time_step = 0.0;
    }
  else
    {
      cumulative_nth_time_step += time_step;
    }

  previous_solve_pressure_velocity_part = solve_pressure_velocity_part;

  std::cout << "   Solving saturation transport equation..." << std::endl;

  assemble_saturation_system ();

  {
    SolverControl solver_control (saturation_matrix.m(),
                                  1e-8*saturation_rhs.l2_norm());
    SolverCG<TrilinosWrappers::Vector> cg (solver_control);

    TrilinosWrappers::PreconditionIC preconditioner;
    preconditioner.initialize (saturation_matrix);

    cg.solve (saturation_matrix, saturation_solution,
              saturation_rhs, preconditioner);


    saturation_constraints.distribute (saturation_solution);

    project_back_saturation ();

    std::cout << "     "
              << solver_control.last_step()
              << " CG iterations for saturation."
              << std::endl;

  }

}

template <int dim>
bool
TwoPhaseFlowProblem<dim>::determine_whether_to_solve_pressure_velocity_part () const
{
  if (timestep_number <= 3)
    return true;

  const QGauss<dim> quadrature_formula(saturation_degree+2);
  const unsigned int   n_q_points = quadrature_formula.size();

  FEValues<dim> fe_values (saturation_fe, quadrature_formula,
                           update_values | update_quadrature_points);

  std::vector<double> old_saturation_after_solving_pressure (n_q_points);
  std::vector<double> present_saturation (n_q_points);

  const RandomMedium::KInverse<dim> k_inverse;
//  const SingleCurvingCrack::KInverse<dim> k_inverse;

  std::vector<Tensor<2,dim> >       k_inverse_values (n_q_points);

  double max_global_aop_indicator = 0.0;

  typename DoFHandler<dim>::active_cell_iterator
    cell = saturation_dof_handler.begin_active(),
    endc = saturation_dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      double max_local_mobility_reciprocal_difference = 0.0;
      double max_local_permeability_inverse_l1_norm = 0.0;

      fe_values.reinit(cell);
      fe_values.get_function_values (nth_saturation_solution_after_solving_pressure_part,
                                     old_saturation_after_solving_pressure);
      fe_values.get_function_values (saturation_solution,
                                     present_saturation);

      k_inverse.value_list (fe_values.get_quadrature_points(),
                            k_inverse_values);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          double mobility_reciprocal_difference = std::fabs( mobility_inverse(present_saturation[q],viscosity)
                                                             -
                                                             mobility_inverse(old_saturation_after_solving_pressure[q],viscosity) );

          max_local_mobility_reciprocal_difference = std::max(max_local_mobility_reciprocal_difference,
                                                              mobility_reciprocal_difference);

          max_local_permeability_inverse_l1_norm = std::max(max_local_permeability_inverse_l1_norm,
                                                            k_inverse_values[q][0][0]);
        }

      max_global_aop_indicator = std::max(max_global_aop_indicator,
                                          (max_local_mobility_reciprocal_difference*max_local_permeability_inverse_l1_norm));
    }

  if ( max_global_aop_indicator > 5.0 )
    {
      return true;
    }
  else
    {
      std::cout << "   Activating adaptive operating splitting" << std::endl;
      return false;
    }
}

template <int dim>
void
TwoPhaseFlowProblem<dim>::
compute_refinement_indicators (Vector<double> &refinement_indicators) const
{

  const QMidpoint<dim> quadrature_formula;
  FEValues<dim> fe_values (saturation_fe, quadrature_formula, update_gradients);
  std::vector<Tensor<1,dim> > grad_saturation (1);

  double max_refinement_indicator = 0.0;

  typename DoFHandler<dim>::active_cell_iterator
    cell = saturation_dof_handler.begin_active(),
    endc = saturation_dof_handler.end();
  for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
    {
      fe_values.reinit(cell);
      fe_values.get_function_grads (predictor_saturation_solution,
                                    grad_saturation);

      refinement_indicators(cell_no)
        = std::log( 1.0 + std::sqrt( grad_saturation[0] *
                                     grad_saturation[0] ) );
      max_refinement_indicator = std::max(max_refinement_indicator,
                                          refinement_indicators(cell_no));
    }

//  std::cout << "max_refinement_indicator =" << max_refinement_indicator << std::endl;
}

template <int dim>
void
TwoPhaseFlowProblem<dim>::
refine_grid (const Vector<double> &refinement_indicators)
{
  const double current_saturation_level = saturation_level +
                                          n_refinement_steps;

  {
    typename DoFHandler<dim>::active_cell_iterator
      cell = saturation_dof_handler.begin_active(),
      endc = saturation_dof_handler.end();

    for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
      {
	cell->clear_coarsen_flag();
	cell->clear_refine_flag();

	if ((cell->level() < current_saturation_level) &&
	    (std::fabs(refinement_indicators(cell_no)) > saturation_value))
	  cell->set_refine_flag();
	else
	    if ((cell->level() > double(n_refinement_steps)) &&
		(std::fabs(refinement_indicators(cell_no)) < 0.75 * saturation_value))
	      cell->set_coarsen_flag();
      }
  }

  triangulation.prepare_coarsening_and_refinement ();

  unsigned int number_of_cells_refine  = 0;
  unsigned int number_of_cells_coarsen = 0;

  {
    typename DoFHandler<dim>::active_cell_iterator
      cell = saturation_dof_handler.begin_active(),
      endc = saturation_dof_handler.end();

    for (; cell!=endc; ++cell)
      if (cell->refine_flag_set())
	    ++number_of_cells_refine;
      else
	if (cell->coarsen_flag_set())
	  ++number_of_cells_coarsen;
  }

  std::cout << "   "
            << number_of_cells_refine
            << " cell(s) are going to be refined."
            << std::endl;
  std::cout << "   "
            << number_of_cells_coarsen
            << " cell(s) are going to be coarsened."
            << std::endl;

  std::cout << std::endl;

  if ( number_of_cells_refine > 0 || number_of_cells_coarsen > 0 )
    {
      std::vector<TrilinosWrappers::Vector> x_saturation (3);
      x_saturation[0] = saturation_solution;
      x_saturation[1] = old_saturation_solution;
      x_saturation[2] = nth_saturation_solution_after_solving_pressure_part;

      std::vector<TrilinosWrappers::BlockVector> x_darcy (2);
      x_darcy[0] = nth_darcy_solution_after_solving_pressure_part;
      x_darcy[1] = n_minus_oneth_darcy_solution_after_solving_pressure_part;

      SolutionTransfer<dim,TrilinosWrappers::Vector> saturation_soltrans(saturation_dof_handler);

      SolutionTransfer<dim,TrilinosWrappers::BlockVector> darcy_soltrans(darcy_dof_handler);


      triangulation.prepare_coarsening_and_refinement();
      saturation_soltrans.prepare_for_coarsening_and_refinement(x_saturation);

      darcy_soltrans.prepare_for_coarsening_and_refinement(x_darcy);

      triangulation.execute_coarsening_and_refinement ();
      setup_dofs ();

      std::vector<TrilinosWrappers::Vector> tmp_saturation (3);
      tmp_saturation[0].reinit (saturation_solution);
      tmp_saturation[1].reinit (saturation_solution);
      tmp_saturation[2].reinit (saturation_solution);
      saturation_soltrans.interpolate(x_saturation, tmp_saturation);

      saturation_solution = tmp_saturation[0];
      old_saturation_solution = tmp_saturation[1];
      nth_saturation_solution_after_solving_pressure_part = tmp_saturation[2];

      std::vector<TrilinosWrappers::BlockVector> tmp_darcy (2);
      tmp_darcy[0].reinit (darcy_solution);
      tmp_darcy[1].reinit (darcy_solution);
      darcy_soltrans.interpolate(x_darcy, tmp_darcy);

      nth_darcy_solution_after_solving_pressure_part = tmp_darcy[0];
      n_minus_oneth_darcy_solution_after_solving_pressure_part = tmp_darcy[1];

      rebuild_saturation_matrix    = true;
    }
  else
    {
      rebuild_saturation_matrix    = false;

      std::vector<unsigned int> darcy_block_component (dim+1,0);
      darcy_block_component[dim] = 1;

      std::vector<unsigned int> darcy_dofs_per_block (2);
      DoFTools::count_dofs_per_block (darcy_dof_handler, darcy_dofs_per_block, darcy_block_component);
      const unsigned int n_u = darcy_dofs_per_block[0],
                         n_p = darcy_dofs_per_block[1],
                         n_s = saturation_dof_handler.n_dofs();

      std::cout << "Number of active cells: "
                << triangulation.n_active_cells()
                << " (on "
                << triangulation.n_levels()
                << " levels)"
                << std::endl
                << "Number of degrees of freedom: "
                << n_u + n_p + n_s
                << " (" << n_u << '+' << n_p << '+'<< n_s <<')'
                << std::endl
                << std::endl;
    }

}


template <int dim>
void TwoPhaseFlowProblem<dim>::output_results ()  const
{
  if ( solve_pressure_velocity_part == false )
    return;

  const FESystem<dim> joint_fe (darcy_fe, 1,
                                saturation_fe, 1);
  DoFHandler<dim> joint_dof_handler (triangulation);
  joint_dof_handler.distribute_dofs (joint_fe);
  Assert (joint_dof_handler.n_dofs() ==
          darcy_dof_handler.n_dofs() + saturation_dof_handler.n_dofs(),
          ExcInternalError());

  Vector<double> joint_solution (joint_dof_handler.n_dofs());

  {
    std::vector<unsigned int> local_joint_dof_indices (joint_fe.dofs_per_cell);
    std::vector<unsigned int> local_darcy_dof_indices (darcy_fe.dofs_per_cell);
    std::vector<unsigned int> local_saturation_dof_indices (saturation_fe.dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
      joint_cell      = joint_dof_handler.begin_active(),
      joint_endc      = joint_dof_handler.end(),
      darcy_cell      = darcy_dof_handler.begin_active(),
      saturation_cell = saturation_dof_handler.begin_active();

    for (; joint_cell!=joint_endc; ++joint_cell, ++darcy_cell, ++saturation_cell)
      {
        joint_cell->get_dof_indices (local_joint_dof_indices);
        darcy_cell->get_dof_indices (local_darcy_dof_indices);
        saturation_cell->get_dof_indices (local_saturation_dof_indices);

        for (unsigned int i=0; i<joint_fe.dofs_per_cell; ++i)
          if (joint_fe.system_to_base_index(i).first.first == 0)
            {
              Assert (joint_fe.system_to_base_index(i).second
                      <
                      local_darcy_dof_indices.size(),
                      ExcInternalError());
              joint_solution(local_joint_dof_indices[i])
                = darcy_solution(local_darcy_dof_indices[joint_fe.system_to_base_index(i).second]);
            }
          else
            {
              Assert (joint_fe.system_to_base_index(i).first.first == 1,
                      ExcInternalError());
              Assert (joint_fe.system_to_base_index(i).second
                      <
                      local_darcy_dof_indices.size(),
                      ExcInternalError());
              joint_solution(local_joint_dof_indices[i])
                = saturation_solution(local_saturation_dof_indices[joint_fe.system_to_base_index(i).second]);
            }

      }
  }
  std::vector<std::string> joint_solution_names;
  switch (dim)
    {
      case 2:
            joint_solution_names.push_back ("u");
            joint_solution_names.push_back ("v");
            break;

      case 3:
            joint_solution_names.push_back ("u");
            joint_solution_names.push_back ("v");
            joint_solution_names.push_back ("w");
            break;

      default:
            Assert (false, ExcNotImplemented());
    }
  joint_solution_names.push_back ("pressure");
  joint_solution_names.push_back ("saturation");

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation
    .push_back (DataComponentInterpretation::component_is_scalar);
  data_component_interpretation
    .push_back (DataComponentInterpretation::component_is_scalar);

  DataOut<dim> data_out;

  data_out.attach_dof_handler (joint_dof_handler);
  data_out.add_data_vector (joint_solution, joint_solution_names,
                            DataOut<dim>::type_dof_data,
                            data_component_interpretation);

  data_out.build_patches ();

  std::string filename = "solution-" +
                         Utilities::int_to_string (timestep_number, 5) + ".tec";
  std::ofstream output (filename.c_str());
  data_out.write_tecplot (output);
}


template <int dim>
void
TwoPhaseFlowProblem<dim>::project_back_saturation ()
{
  for (unsigned int i=0; i<saturation_solution.size(); ++i)
    if (saturation_solution(i) < 0)
      saturation_solution(i) = 0;
    else
      if (saturation_solution(i) > 1)
        saturation_solution(i) = 1;
}


template <int dim>
double
TwoPhaseFlowProblem<dim>::get_maximal_velocity () const
{
  QGauss<dim>   quadrature_formula(darcy_degree+2);
  const unsigned int   n_q_points
    = quadrature_formula.size();

  FEValues<dim> darcy_fe_values (darcy_fe, quadrature_formula,
                                 update_values);
  std::vector<Vector<double> > darcy_solution_values(n_q_points,
                                                     Vector<double>(dim+1));
  double max_velocity = 0;

  typename DoFHandler<dim>::active_cell_iterator
    cell = darcy_dof_handler.begin_active(),
    endc = darcy_dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      darcy_fe_values.reinit (cell);
      darcy_fe_values.get_function_values (darcy_solution, darcy_solution_values);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          Tensor<1,dim> velocity;
          for (unsigned int i=0; i<dim; ++i)
            velocity[i] = darcy_solution_values[q](i);

          max_velocity = std::max (max_velocity,
                                   velocity.norm());
        }
    }

  return max_velocity;
}


template <int dim>
std::pair<double,double>
TwoPhaseFlowProblem<dim>::get_extrapolated_saturation_range () const
{
  const QGauss<dim>  quadrature_formula(saturation_degree+2);
  const unsigned int n_q_points = quadrature_formula.size();

  FEValues<dim> fe_values (saturation_fe, quadrature_formula,
                           update_values);
  std::vector<double> old_saturation_values(n_q_points);
  std::vector<double> old_old_saturation_values(n_q_points);

  if (timestep_number != 0)
    {
      double min_saturation = (1. + time_step/old_time_step) *
                              old_saturation_solution.linfty_norm()
                              +
                              time_step/old_time_step *
                              old_old_saturation_solution.linfty_norm(),
             max_saturation = -min_saturation;

      typename DoFHandler<dim>::active_cell_iterator
        cell = saturation_dof_handler.begin_active(),
        endc = saturation_dof_handler.end();
      for (; cell!=endc; ++cell)
        {
          fe_values.reinit (cell);
          fe_values.get_function_values (old_saturation_solution,
                                         old_saturation_values);
          fe_values.get_function_values (old_old_saturation_solution,
                                         old_old_saturation_values);

          for (unsigned int q=0; q<n_q_points; ++q)
            {
              const double saturation =
                (1. + time_step/old_time_step) * old_saturation_values[q]-
                time_step/old_time_step * old_old_saturation_values[q];

              min_saturation = std::min (min_saturation, saturation);
              max_saturation = std::max (max_saturation, saturation);
            }
        }

      return std::make_pair(min_saturation, max_saturation);
    }
  else
    {
      double min_saturation = old_saturation_solution.linfty_norm(),
             max_saturation = -min_saturation;

      typename DoFHandler<dim>::active_cell_iterator
        cell = saturation_dof_handler.begin_active(),
        endc = saturation_dof_handler.end();
      for (; cell!=endc; ++cell)
        {
          fe_values.reinit (cell);
          fe_values.get_function_values (old_saturation_solution,
                                         old_saturation_values);

          for (unsigned int q=0; q<n_q_points; ++q)
            {
              const double saturation = old_saturation_values[q];

              min_saturation = std::min (min_saturation, saturation);
              max_saturation = std::max (max_saturation, saturation);
            }
        }

      return std::make_pair(min_saturation, max_saturation);
    }
}

template <int dim>
double
TwoPhaseFlowProblem<dim>::
compute_viscosity (const std::vector<double>          &old_saturation,
                   const std::vector<double>          &old_old_saturation,
                   const std::vector<Tensor<1,dim> >  &old_saturation_grads,
                   const std::vector<Tensor<1,dim> >  &old_old_saturation_grads,
                   const std::vector<Vector<double> > &present_darcy_values,
                   const double                        global_u_infty,
                   const double                        global_S_variation,
                   const double                        global_Omega_diameter,
                   const double                        cell_diameter,
                   const double                        old_time_step,
                   const double                        viscosity)
{
  const double beta = 0.08 * dim;
  const double alpha = 1;

  if (global_u_infty == 0)
    return 5e-3 * cell_diameter;

  const unsigned int n_q_points = old_saturation.size();

  double max_residual = 0;
  double max_velocity = 0;

  for (unsigned int q=0; q < n_q_points; ++q)
    {
      Tensor<1,dim> u;
      for (unsigned int d=0; d<dim; ++d)
        u[d] = present_darcy_values[q](d);

      const double dS_dt = (old_saturation[q] - old_old_saturation[q])
                           / old_time_step;

      const double dF_dS = get_fractional_flow_derivative ((old_saturation[q] + old_old_saturation[q]) / 2.0,
                                                           viscosity);

      const double u_grad_S = u * dF_dS *
                              (old_saturation_grads[q] + old_old_saturation_grads[q]) / 2.0;

      const double residual
        = std::abs((dS_dt + u_grad_S) *
                   std::pow((old_saturation[q]+old_old_saturation[q]) / 2,
                            alpha-1.));

      max_residual = std::max (residual,        max_residual);
      max_velocity = std::max (std::sqrt (u*u), max_velocity);
    }

  const double global_scaling = global_u_infty * global_S_variation /
                                std::pow(global_Omega_diameter, alpha - 2.);

  return (beta *
          max_velocity *
          std::min (cell_diameter,
                    std::pow(cell_diameter,alpha) *
                    max_residual / global_scaling));
}


template <int dim>
void TwoPhaseFlowProblem<dim>::run ()
{
  unsigned int pre_refinement_step = 0;

  GridGenerator::hyper_cube (triangulation, 0, 1);
  triangulation.refine_global (n_refinement_steps);

  setup_dofs ();

  start_time_iteration:

  VectorTools::project (saturation_dof_handler,
                        saturation_constraints,
                        QGauss<dim>(saturation_degree+2),
                        SaturationInitialValues<dim>(),
                        old_saturation_solution);

  timestep_number = 0;
  double time = 0;

  do
    {
      std::cout << "Timestep " << timestep_number
                << ":  t=" << time
                << ", dt=" << time_step
                << std::endl;

      solve ();

      output_results ();

      solve_pressure_velocity_part = false;

      if ((timestep_number == 0) &&
          (pre_refinement_step < saturation_level))
        {
          predictor_saturation_solution = saturation_solution;
          predictor_saturation_solution.sadd (2.0, -1.0, old_saturation_solution);
          Vector<double> refinement_indicators (triangulation.n_active_cells());
          compute_refinement_indicators(refinement_indicators);
          refine_grid(refinement_indicators);
          ++pre_refinement_step;
          goto start_time_iteration;
        }
      else
        {
          predictor_saturation_solution = saturation_solution;
          predictor_saturation_solution.sadd (2.0, -1.0, old_saturation_solution);
          Vector<double> refinement_indicators (triangulation.n_active_cells());
          compute_refinement_indicators(refinement_indicators);
          refine_grid(refinement_indicators);
        }

      time += time_step;
      ++timestep_number;

      old_old_saturation_solution = old_saturation_solution;
      old_saturation_solution = saturation_solution;

    }
  while (time <= 250);
}


int main ()
{
  try
    {
      deallog.depth_console (0);

      TwoPhaseFlowProblem<3> two_phase_flow_problem(1);
      two_phase_flow_problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
