/* $Id$ */
/* Author: Wolfgang Bangerth, Texas A&M University, 2007 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2007, 2008 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */



#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function.h>
#include <base/utilities.h>

#include <lac/block_vector.h>
#include <lac/full_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/solver_gmres.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/sparse_direct.h>
#include <lac/sparse_ilu.h>
#include <lac/block_matrix_array.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_tools.h>
#include <grid/grid_refinement.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_constraints.h>

#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>

#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/derivative_approximation.h>
#include <numerics/solution_transfer.h>

#include <fstream>
#include <sstream>

using namespace dealii;


template <int dim>
 struct InnerPreconditioner;

 template <>
 struct InnerPreconditioner<2>
 {
     typedef SparseDirectUMFPACK type;
 };

 template <>
 struct InnerPreconditioner<3>
 {
     typedef SparseILU<double> type;
 };


template <int dim>
class BoussinesqFlowProblem
{
  public:
    BoussinesqFlowProblem (const unsigned int degree);
    void run ();

  private:
    void setup_dofs (const bool setup_matrices);
    void assemble_system ();
    void assemble_rhs_T ();
    double get_maximal_velocity () const;
    void solve ();
    void output_results () const;
    void refine_mesh ();

    const unsigned int   degree;

    Triangulation<dim>   triangulation;
    FESystem<dim>        fe;
    DoFHandler<dim>      dof_handler;

    ConstraintMatrix     hanging_node_constraints;

    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    double time_step;
    unsigned int timestep_number;

    BlockVector<double> solution;
    BlockVector<double> old_solution;
    BlockVector<double> system_rhs;

    boost::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner;
    boost::shared_ptr<SparseILU<double> > Mp_preconditioner;

    bool rebuild_matrices;
    bool rebuild_preconditioner;
};





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
PressureBoundaryValues<dim>::value (const Point<dim>  &/*p*/,
                                    const unsigned int /*component*/) const
{
  return 0;
}



template <int dim>
class TemperatureBoundaryValues : public Function<dim>
{
  public:
    TemperatureBoundaryValues () : Function<dim>(1) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
};



template <int dim>
double
TemperatureBoundaryValues<dim>::value (const Point<dim> &p,
                                      const unsigned int /*component*/) const
{
//TODO: leftover from olden times. replace by something sensible once we have
//diffusion in the temperature field
  if (p[0] == 0)
    return 1;
  else
    return 0;
}




template <int dim>
class InitialValues : public Function<dim>
{
  public:
    InitialValues () : Function<dim>(dim+2) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
};


template <int dim>
double
InitialValues<dim>::value (const Point<dim>  &,
                           const unsigned int) const
{
  return 0;
}


template <int dim>
void
InitialValues<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = InitialValues<dim>::value (p, c);
}



template <int dim>
class RightHandSide : public Function<dim>
{
  public:
    RightHandSide () : Function<dim>(dim+2) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
};


template <int dim>
double
RightHandSide<dim>::value (const Point<dim>  &p,
                           const unsigned int component) const
{
  if (component == dim+1)
    return ((p.distance (Point<dim>(.3,.1)) < 1./32)
	    ||
	    (p.distance (Point<dim>(.45,.1)) < 1./32)
	    ||
	    (p.distance (Point<dim>(.75,.1)) < 1./32)
	    ?
	    1
	    :
	    0);
  else
    return 0;
}


template <int dim>
void
RightHandSide<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = RightHandSide<dim>::value (p, c);
}




template <class Matrix, class Preconditioner>
class InverseMatrix : public Subscriptor
{
  public:
    InverseMatrix (const Matrix         &m,
		   const Preconditioner &preconditioner);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const Matrix> matrix;
    const Preconditioner &preconditioner;
};


template <class Matrix, class Preconditioner>
InverseMatrix<Matrix,Preconditioner>::InverseMatrix (const Matrix &m,
						     const Preconditioner &preconditioner)
                :
                matrix (&m),
		preconditioner (preconditioner)
{}



template <class Matrix, class Preconditioner>
void InverseMatrix<Matrix,Preconditioner>::vmult (Vector<double>       &dst,
						  const Vector<double> &src) const
{
  SolverControl solver_control (src.size(), 1e-6*src.l2_norm());
  SolverCG<> cg (solver_control);

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


				// This is the implementation
				// of the Schur complement
				// preconditioner as described
				// in the section on improved
				// solvers in step-22. See there
				// for more explanation of the
				// method.
template <class PreconditionerA, class PreconditionerMp>
class BlockSchurPreconditioner : public Subscriptor
{
  public:
    BlockSchurPreconditioner (const BlockSparseMatrix<double>         &S,
          const InverseMatrix<SparseMatrix<double>,PreconditionerMp>  &Mpinv,
          const PreconditionerA &Apreconditioner);

  void vmult (BlockVector<double>       &dst,
              const BlockVector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    const SmartPointer<const InverseMatrix<SparseMatrix<double>,
                       PreconditionerMp > > m_inverse;
    const PreconditionerA &a_preconditioner;

    mutable Vector<double> tmp;

};

template <class PreconditionerA, class PreconditionerMp>
BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::BlockSchurPreconditioner(
          const BlockSparseMatrix<double>                            &S,
          const InverseMatrix<SparseMatrix<double>,PreconditionerMp> &Mpinv,
          const PreconditionerA &Apreconditioner
          )
                :
                system_matrix           (&S),
                m_inverse               (&Mpinv),
                a_preconditioner        (Apreconditioner),
                tmp                     (S.block(1,1).m())
{
}

template <class PreconditionerA, class PreconditionerMp>
void BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::vmult (
                                     BlockVector<double>       &dst,
                                     const BlockVector<double> &src) const
{
  a_preconditioner.vmult (dst.block(0), src.block(0));
  system_matrix->block(1,0).residual(tmp, dst.block(0), src.block(1));
  tmp *= -1;
  m_inverse->vmult (dst.block(1), tmp);
}



template <int dim>
BoussinesqFlowProblem<dim>::BoussinesqFlowProblem (const unsigned int degree)
                :
                degree (degree),
                fe (FE_Q<dim>(degree+1), dim,
                    FE_Q<dim>(degree), 1,
                    FE_DGQ<dim>(degree-1), 1),
                dof_handler (triangulation),
                time_step (0),
		rebuild_matrices (true),
		rebuild_preconditioner (true)
{}




template <int dim>
void BoussinesqFlowProblem<dim>::setup_dofs (const bool setup_matrices)
{
  dof_handler.distribute_dofs (fe);
  DoFRenumbering::component_wise (dof_handler);

  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
  std::set<unsigned char> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert (0);
  VectorTools::compute_no_normal_flux_constraints (dof_handler, 0,
						   no_normal_flux_boundaries,
						   hanging_node_constraints);
  hanging_node_constraints.close ();

  std::vector<unsigned int> dofs_per_component (dim+2);
  DoFTools::count_dofs_per_component (dof_handler, dofs_per_component);
  const unsigned int n_u = dofs_per_component[0] * dim,
                     n_p = dofs_per_component[dim],
                     n_T = dofs_per_component[dim+1];

  std::cout << "Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl
            << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p << '+'<< n_T <<')'
            << std::endl
            << std::endl;

  const unsigned int
    n_couplings = dof_handler.max_couplings_between_dofs();

  if (setup_matrices == true)
    {
      system_matrix.clear ();

      sparsity_pattern.reinit (3,3);
      sparsity_pattern.block(0,0).reinit (n_u, n_u, n_couplings);
      sparsity_pattern.block(1,0).reinit (n_p, n_u, n_couplings);
      sparsity_pattern.block(2,0).reinit (n_T, n_u, n_couplings);
      sparsity_pattern.block(0,1).reinit (n_u, n_p, n_couplings);
      sparsity_pattern.block(1,1).reinit (n_p, n_p, n_couplings);
      sparsity_pattern.block(2,1).reinit (n_T, n_p, n_couplings);
      sparsity_pattern.block(0,2).reinit (n_u, n_T, n_couplings);
      sparsity_pattern.block(1,2).reinit (n_p, n_T, n_couplings);
      sparsity_pattern.block(2,2).reinit (n_T, n_T, n_couplings);

      sparsity_pattern.collect_sizes();


      DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
      hanging_node_constraints.condense (sparsity_pattern);
      sparsity_pattern.compress();

      system_matrix.reinit (sparsity_pattern);
    }

  solution.reinit (3);
  solution.block(0).reinit (n_u);
  solution.block(1).reinit (n_p);
  solution.block(2).reinit (n_T);
  solution.collect_sizes ();

  old_solution.reinit (3);
  old_solution.block(0).reinit (n_u);
  old_solution.block(1).reinit (n_p);
  old_solution.block(2).reinit (n_T);
  old_solution.collect_sizes ();

  system_rhs.reinit (3);
  system_rhs.block(0).reinit (n_u);
  system_rhs.block(1).reinit (n_p);
  system_rhs.block(2).reinit (n_T);
  system_rhs.collect_sizes ();
}



template <int dim>
void BoussinesqFlowProblem<dim>::assemble_system ()
{
  if (rebuild_matrices == true)
    system_matrix=0;

  system_rhs=0;

  QGauss<dim>   quadrature_formula(degree+2);
  QGauss<dim-1> face_quadrature_formula(degree+2);

  FEValues<dim> fe_values (fe, quadrature_formula,
			   update_values    |
			   update_quadrature_points  |
			   update_JxW_values |
			   (rebuild_matrices == true
			    ?
			    update_gradients
			    :
			    UpdateFlags(0)));
  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
                                    update_values    | update_normal_vectors |
                                    update_quadrature_points  | update_JxW_values);

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;

  const unsigned int   n_q_points      = quadrature_formula.size();
  const unsigned int   n_face_q_points = face_quadrature_formula.size();

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  const PressureBoundaryValues<dim> pressure_boundary_values;

  std::vector<double>               boundary_values (n_face_q_points);

  std::vector<Vector<double> >      old_solution_values(n_q_points, Vector<double>(dim+2));
  std::vector<std::vector<Tensor<1,dim> > >  old_solution_grads(n_q_points,
                                                                std::vector<Tensor<1,dim> > (dim+2));

  const double Raleigh_number = 10;

  std::vector<Tensor<1,dim> >          phi_u       (dofs_per_cell);
  std::vector<SymmetricTensor<2,dim> > phi_grads_u (dofs_per_cell);
  std::vector<double>                  div_phi_u   (dofs_per_cell);
  std::vector<double>                  phi_p       (dofs_per_cell);
  std::vector<double>                  phi_T       (dofs_per_cell);
  std::vector<Tensor<1,dim> >          grad_phi_T  (dofs_per_cell);

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);
  const FEValuesExtractors::Scalar temperature (dim+1);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;

      fe_values.get_function_values (old_solution, old_solution_values);

      for (unsigned int q=0; q<n_q_points; ++q)
	{
	  const double old_temperature = old_solution_values[q](dim+1);

				// Extract the basis relevant
				// terms in the inner products
				// once in advance as shown
				// in step-22. This accelerates
				// the assembly process,
	  for (unsigned int k=0; k<dofs_per_cell; ++k)
	    {
	      phi_u[k] = fe_values[velocities].value (k,q);
	      if (rebuild_matrices)
	        {
		  phi_grads_u[k] = fe_values[velocities].symmetric_gradient(k,q);
		  div_phi_u[k]   = fe_values[velocities].divergence (k, q);
		  phi_p[k]       = fe_values[pressure].value (k, q);
		  phi_T[k]       = fe_values[temperature].value (k, q);
		  grad_phi_T[k]  = fe_values[temperature].gradient (k, q);
		}
	    }

	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {

	      const Tensor<1,dim> phi_i_u = fe_values[velocities].value (i, q);

	      if (rebuild_matrices)
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  local_matrix(i,j) += (phi_grads_u[i] * phi_grads_u[j]
					- div_phi_u[i] * phi_p[j]
					- phi_p[i] * div_phi_u[j]
					+ phi_p[i] * phi_p[j]
					+ phi_T[i] * phi_T[j])
				       * fe_values.JxW(q);

	      const Point<dim> gravity (0,1);

	      local_rhs(i) += (Raleigh_number *
			       gravity * phi_u[i] * old_temperature)*
			      fe_values.JxW(q);
          }
	}


				// The assembly of the face
				// cells which enters the
				// right hand sides cannot
				// be accelerated with the
				// above technique, since
				// all the basis functions are
				// only evaluated once.
      for (unsigned int face_no=0;
           face_no<GeometryInfo<dim>::faces_per_cell;
           ++face_no)
        if (cell->at_boundary(face_no))
          {
            fe_face_values.reinit (cell, face_no);

            pressure_boundary_values
              .value_list (fe_face_values.get_quadrature_points(),
                           boundary_values);

            for (unsigned int q=0; q<n_face_q_points; ++q)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  const Tensor<1,dim>
                    phi_i_u = fe_face_values[velocities].value (i, q);

                  local_rhs(i) += -(phi_i_u *
                                    fe_face_values.normal_vector(q) *
                                    boundary_values[q] *
                                    fe_face_values.JxW(q));
                }
          }

      cell->get_dof_indices (local_dof_indices);

      if (rebuild_matrices == true)
	{
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      system_matrix.add (local_dof_indices[i],
				 local_dof_indices[j],
				 local_matrix(i,j));
	}

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += local_rhs(i);
    }

  if (rebuild_matrices == true)
    hanging_node_constraints.condense (system_matrix);

  hanging_node_constraints.condense (system_rhs);

  if (rebuild_matrices == true)
    {
//       std::map<unsigned int,double> boundary_values;

//       typename DoFHandler<dim>::active_cell_iterator
// 	cell = dof_handler.begin_active(),
// 	emdc = dof_handler.end();
//       for (; cell!=endc; ++cell)
// 	for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
// 	  if (cell->vertex(v).distance(dim == 2
// 				       ?
// 				       Point<dim>(0,-1)
// 				       :
// 				       Point<dim>(0,0,-1)) < 1e-6)
// 	    {
// 	      std::cout << "Found cell and vertex: " << cell << ' '
// 			<< v << std::endl;

// 	      boundary_values[cell->vertex_dof_index(v,0)] = 0;
// 	      break;
// 	    }

//      std::vector<bool> component_mask (dim+2, true);
//       component_mask[dim] = component_mask[dim+1] = false;
//       VectorTools::interpolate_boundary_values (dof_handler,
// 						0,
// 						ZeroFunction<dim>(dim+2),
// 						boundary_values,
// 						component_mask);

//       MatrixTools::apply_boundary_values (boundary_values,
// 					  system_matrix,
// 					  solution,
// 					  system_rhs);
    }

  if (rebuild_preconditioner == true)
    {
      Assert (rebuild_matrices == true,
	      ExcMessage ("There is no point in rebuilding the preconditioner "
			  "without a rebuilt matrix!"));

      std::cout << "   Rebuilding preconditioner..." << std::flush;

				// Rebuild the preconditioner
				// for the velocity-velocity
				// block (0,0)
      A_preconditioner
	= boost::shared_ptr<typename InnerPreconditioner<dim>::type>
		(new typename InnerPreconditioner<dim>::type());
      A_preconditioner->initialize (system_matrix.block(0,0),
		typename InnerPreconditioner<dim>::type::AdditionalData());

				// Rebuild the preconditioner
				// for the pressure-pressure
				// block (1,1)
      Mp_preconditioner
	= boost::shared_ptr<SparseILU<double> >
		(new SparseILU<double>);
      Mp_preconditioner->initialize (system_matrix.block(1,1),
				     SparseILU<double>::AdditionalData());

      std::cout << std::endl;

      rebuild_preconditioner = false;
    }

  rebuild_matrices = false;
}






template <int dim>
void BoussinesqFlowProblem<dim>::assemble_rhs_T ()
{
  QGauss<dim>   quadrature_formula(degree+2);
  QGauss<dim-1> face_quadrature_formula(degree+2);
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    | update_gradients |
                           update_quadrature_points  | update_JxW_values);
  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
                                    update_values    | update_normal_vectors |
                                    update_quadrature_points  | update_JxW_values);
  FESubfaceValues<dim> fe_subface_values (fe, face_quadrature_formula,
					  update_values    | update_normal_vectors |
					  update_JxW_values);
  FEFaceValues<dim> fe_face_values_neighbor (fe, face_quadrature_formula,
                                             update_values);
  FESubfaceValues<dim> fe_subface_values_neighbor (fe, face_quadrature_formula,
						   update_values);

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  const unsigned int   n_q_points      = quadrature_formula.size();
  const unsigned int   n_face_q_points = face_quadrature_formula.size();

  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<Vector<double> > old_solution_values(n_q_points, Vector<double>(dim+2));

  std::vector<Vector<double> > old_solution_values_face(n_face_q_points, Vector<double>(dim+2));
  std::vector<Vector<double> > old_solution_values_face_neighbor(n_face_q_points, Vector<double>(dim+2));
  std::vector<Vector<double> > present_solution_values(n_q_points, Vector<double>(dim+2));
  std::vector<Vector<double> > present_solution_values_face(n_face_q_points, Vector<double>(dim+2));

  std::vector<std::vector<Tensor<1,dim> > >
    present_solution_grads(n_q_points,
			   std::vector<Tensor<1,dim> >(dim+2));


  std::vector<double> neighbor_temperature (n_face_q_points);
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  TemperatureBoundaryValues<dim> temperature_boundary_values;
  const FEValuesExtractors::Scalar temperature (dim+1);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  for (; cell!=endc; ++cell)
    {
      local_rhs = 0;
      fe_values.reinit (cell);

      fe_values.get_function_values (old_solution, old_solution_values);
      fe_values.get_function_values (solution, present_solution_values);
      fe_values.get_function_gradients (solution, present_solution_grads);

      for (unsigned int q=0; q<n_q_points; ++q)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const double old_T = old_solution_values[q](dim+1);
            Tensor<1,dim> present_u;
            for (unsigned int d=0; d<dim; ++d)
              present_u[d] = present_solution_values[q](d);

	    double present_div_u = 0;
            for (unsigned int d=0; d<dim; ++d)
              present_div_u += present_solution_grads[q][d][d];

            const double        phi_i_T      = fe_values[temperature].value (i, q);
            const Tensor<1,dim> grad_phi_i_T = fe_values[temperature].gradient (i, q);

	    const Point<dim> p = fe_values.quadrature_point(q);

            local_rhs(i) += (time_step *
                             old_T *
                             (present_u *
			      grad_phi_i_T
			      +
			      present_div_u *
			      phi_i_T)
                             +
                             old_T * phi_i_T
			     +
			     time_step *
			     RightHandSide<dim>().value (p, dim+1)
			     * phi_i_T)
                            *
                            fe_values.JxW(q);
          }


//TODO: unify the code that actually does the assembly down below
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
           ++face_no)
        if (cell->at_boundary(face_no)
	    ||
	    ((cell->neighbor(face_no)->has_children() == false)
	     &&
	     (cell->neighbor(face_no)->level() == cell->level())))
	  {
					     // cell either at
					     // boundary or with a
					     // neighbor that has the
					     // same refinement level
					     // and is not further
					     // refined
	    fe_face_values.reinit (cell, face_no);

	    fe_face_values.get_function_values (old_solution,
						old_solution_values_face);
	    fe_face_values.get_function_values (solution,
						present_solution_values_face);

	    if (cell->at_boundary(face_no))
	      temperature_boundary_values
		.value_list (fe_face_values.get_quadrature_points(),
			     neighbor_temperature);
	    else
	      {
		const typename DoFHandler<dim>::active_cell_iterator
		  neighbor = cell->neighbor(face_no);

		fe_face_values_neighbor.reinit (neighbor,
						cell->neighbor_of_neighbor(face_no));

		fe_face_values_neighbor
		  .get_function_values (old_solution,
					old_solution_values_face_neighbor);

		for (unsigned int q=0; q<n_face_q_points; ++q)
		  neighbor_temperature[q] = old_solution_values_face_neighbor[q](dim+1);
	      }

	    for (unsigned int q=0; q<n_face_q_points; ++q)
	      {
		Tensor<1,dim> present_u_face;
		for (unsigned int d=0; d<dim; ++d)
		  present_u_face[d] = present_solution_values_face[q](d);

		const double normal_flux = present_u_face *
					   fe_face_values.normal_vector(q);

		const bool is_outflow_q_point = (normal_flux >= 0);

		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  local_rhs(i) -= time_step *
				  normal_flux *
				  (is_outflow_q_point == true
				   ?
				   old_solution_values_face[q](dim+1)
				   :
				   neighbor_temperature[q]) *
				  fe_face_values[temperature].value (i,q) *
				  fe_face_values.JxW(q);
	      }
	  }
	else
	  if (cell->neighbor(face_no)->has_children())
	    {
					       // neighbor is further
					       // refined. loop over
					       // all sub faces
	      for (unsigned int subface_no=0;
		   subface_no<GeometryInfo<dim>::subfaces_per_face;
		   ++subface_no)
		{
		  fe_subface_values.reinit (cell, face_no, subface_no);

		  fe_subface_values.get_function_values (old_solution,
							 old_solution_values_face);
		  fe_subface_values.get_function_values (solution,
							 present_solution_values_face);

		  const typename DoFHandler<dim>::active_cell_iterator
		    neighbor = cell->neighbor_child_on_subface (face_no, subface_no);

		  fe_face_values_neighbor.reinit (neighbor,
						  cell->neighbor_of_neighbor(face_no));

		  fe_face_values_neighbor
		    .get_function_values (old_solution,
					  old_solution_values_face_neighbor);

		  for (unsigned int q=0; q<n_face_q_points; ++q)
		    neighbor_temperature[q] = old_solution_values_face_neighbor[q](dim+1);

		  for (unsigned int q=0; q<n_face_q_points; ++q)
		    {
		      Tensor<1,dim> present_u_face;
		      for (unsigned int d=0; d<dim; ++d)
			present_u_face[d] = present_solution_values_face[q](d);

		      const double normal_flux = present_u_face *
						 fe_subface_values.normal_vector(q);

		      const bool is_outflow_q_point = (normal_flux >= 0);

		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			local_rhs(i) -= time_step *
					normal_flux *
					(is_outflow_q_point == true
					 ?
					 old_solution_values_face[q](dim+1)
					 :
					 neighbor_temperature[q]) *
					fe_face_values[temperature].value (i,q) *
					fe_face_values.JxW(q);
		    }
		}
	    }
	  else
	    {
					       // neighbor is less
					       // refined. we need to
					       // use a subface values
					       // object for the
					       // neighbor's subface
	      fe_face_values.reinit (cell, face_no);

	      fe_face_values.get_function_values (old_solution, old_solution_values_face);
	      fe_face_values.get_function_values (solution, present_solution_values_face);

	      const typename DoFHandler<dim>::active_cell_iterator
		neighbor = cell->neighbor (face_no);

	      const std::pair<unsigned int, unsigned int> faceno_subfaceno=
		cell->neighbor_of_coarser_neighbor(face_no);
	      const unsigned int neighbor_face_no    = faceno_subfaceno.first,
				 neighbor_subface_no = faceno_subfaceno.second;

	      Assert (neighbor->neighbor_child_on_subface (neighbor_face_no,
							   neighbor_subface_no)
		      == cell,
		      ExcInternalError());

	      fe_subface_values_neighbor.reinit (neighbor,
						 neighbor_face_no,
						 neighbor_subface_no);

	      fe_subface_values_neighbor
		.get_function_values (old_solution,
				      old_solution_values_face_neighbor);

	      for (unsigned int q=0; q<n_face_q_points; ++q)
		neighbor_temperature[q] = old_solution_values_face_neighbor[q](dim+1);

	      for (unsigned int q=0; q<n_face_q_points; ++q)
		{
		  Tensor<1,dim> present_u_face;
		  for (unsigned int d=0; d<dim; ++d)
		    present_u_face[d] = present_solution_values_face[q](d);

		  const double normal_flux = present_u_face *
					     fe_face_values.normal_vector(q);

		  const bool is_outflow_q_point = (normal_flux >= 0);

		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    local_rhs(i) -= time_step *
				    normal_flux *
				    (is_outflow_q_point == true
				     ?
				     old_solution_values_face[q](dim+1)
				     :
				     neighbor_temperature[q]) *
				    fe_face_values[temperature].value (i,q) *
				    fe_face_values.JxW(q);
		}
	    }

      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += local_rhs(i);
    }
}




template <int dim>
void BoussinesqFlowProblem<dim>::solve ()
{
  solution = old_solution;

				// Use the BlockMatrixArray structure
				// for extracting only the upper left
				// 2x2 blocks from the matrix that will
				// be used for the solution of the
				// blocked system.
  {
    GrowingVectorMemory<Vector<double> > simple_mem;
    BlockMatrixArray<double> stokes_submatrix(2, 2, simple_mem);

    stokes_submatrix.enter(system_matrix.block(0,0),0,0);
    stokes_submatrix.enter(system_matrix.block(0,1),0,1);
    stokes_submatrix.enter(system_matrix.block(1,0),1,0);

				// Define some temporary vectors
				// for the solution process.
				// TODO: Can we somhow avoid copying
				// the vectors back and forth? I.e.
				// accessing the block vectors in a
				// similar way as the matrix with the
				// BlockMatrixArray class?
    std::vector<unsigned int> block_sizes(2);
    block_sizes[0] = solution.block(0).size();
    block_sizes[1] = solution.block(1).size();

    BlockVector<double> up_rhs(block_sizes);
    BlockVector<double> up(block_sizes);

    up_rhs.block(0) = system_rhs.block(0);
    up_rhs.block(1) = system_rhs.block(1);

				// Set up inverse matrix for
				// pressure mass matrix
    InverseMatrix<SparseMatrix<double>,SparseILU<double> >
      mp_inverse (system_matrix.block(1,1), *Mp_preconditioner);

				// Set up block Schur preconditioner
    BlockSchurPreconditioner<typename InnerPreconditioner<dim>::type,
                           SparseILU<double> >
      preconditioner (system_matrix, mp_inverse, *A_preconditioner);

				// Set up GMRES solver and
				// solve.
    SolverControl solver_control (system_matrix.m(),
                                    1e-6*system_rhs.l2_norm());

    SolverGMRES<BlockVector<double> > gmres(solver_control,
			SolverGMRES<BlockVector<double> >::AdditionalData(100));

    gmres.solve(stokes_submatrix, up, up_rhs,
                preconditioner);

				// Produce a constistent solution field
    hanging_node_constraints.distribute (up);

    solution.block(0) = up.block(0);
    solution.block(1) = up.block(1);
  }
				   // for DGQ1 needs to be /15
  time_step = GridTools::minimal_cell_diameter(triangulation) /
              std::max (get_maximal_velocity(), .05) / 2;

  assemble_rhs_T ();
  {

    SolverControl solver_control (system_matrix.block(2,2).m(),
                                  1e-8*system_rhs.block(2).l2_norm());
    SolverCG<>   cg (solver_control);
    PreconditionJacobi<> preconditioner;
    preconditioner.initialize (system_matrix.block(2,2));

    try
      {
	cg.solve (system_matrix.block(2,2), solution.block(2), system_rhs.block(2),
		  preconditioner);
      }
    catch (...)
      {
	abort ();
      }

				     // produce a consistent temperature field
    hanging_node_constraints.distribute (solution);

    std::cout << "   "
              << solver_control.last_step()
              << " CG iterations for temperature."
              << std::endl;
    std::cout << "   Max temperature: "
	      << *std::max_element (solution.block(2).begin(),
				    solution.block(2).end())
	      << std::endl;
  }
}



template <int dim>
void BoussinesqFlowProblem<dim>::output_results ()  const
{
  if (timestep_number % 10 != 0)
    return;

  std::vector<std::string> solution_names (dim, "velocity");
  solution_names.push_back ("p");
  solution_names.push_back ("T");

  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim+2, DataComponentInterpretation::component_is_scalar);
  for (unsigned int i=0; i<dim; ++i)
    data_component_interpretation[i]
      = DataComponentInterpretation::component_is_part_of_vector;

  data_out.add_data_vector (solution, solution_names,
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);

  data_out.build_patches ();

  std::ostringstream filename;
  filename << "solution-" << Utilities::int_to_string(timestep_number, 4) << ".vtk";

  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);
}



template <int dim>
void
BoussinesqFlowProblem<dim>::refine_mesh ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

//TODO do this better
  DerivativeApproximation::approximate_gradient (dof_handler,
						 old_solution,
						 estimated_error_per_cell,
						 dim+1);

  typename Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
  for (unsigned int cell_index=0; cell!=endc; ++cell, ++cell_index)
    estimated_error_per_cell(cell_index) *= cell->diameter();

  GridRefinement::refine_and_coarsen_fixed_fraction (triangulation,
						     estimated_error_per_cell,
						     0.3, 0.03,
						     static_cast<unsigned int>
						     (triangulation.n_active_cells()*1.1));

  SolutionTransfer<dim, double> soltrans(dof_handler);

  triangulation.prepare_coarsening_and_refinement();

  Vector<double> x_old_solution (dof_handler.n_dofs());
  x_old_solution = old_solution;

  soltrans.prepare_for_coarsening_and_refinement(x_old_solution);

  triangulation.execute_coarsening_and_refinement ();
  setup_dofs (true);

  Vector<double> tmp (dof_handler.n_dofs());
  soltrans.interpolate(x_old_solution, tmp);

  rebuild_matrices       = true;
  rebuild_preconditioner = true;

  old_solution = tmp;
}



template <int dim>
double
BoussinesqFlowProblem<dim>::get_maximal_velocity () const
{
  QGauss<dim>   quadrature_formula(degree+2);
  const unsigned int   n_q_points
    = quadrature_formula.size();

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values);
  std::vector<Vector<double> > solution_values(n_q_points,
                                               Vector<double>(dim+2));
  double max_velocity = 0;

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values (solution, solution_values);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          Tensor<1,dim> velocity;
          for (unsigned int i=0; i<dim; ++i)
            velocity[i] = solution_values[q](i);

          max_velocity = std::max (max_velocity,
                                   velocity.norm());
        }
    }

  return max_velocity;
}



template <int dim>
void BoussinesqFlowProblem<dim>::run ()
{
  switch (dim)
    {
      case 2:
      {
// 	GridGenerator::hyper_ball (triangulation);

// 	static const HyperBallBoundary<dim> boundary;
// 	triangulation.set_boundary (0, boundary);

	GridGenerator::hyper_cube (triangulation);

	triangulation.refine_global (6);

	break;
      }

      case 3:
      {
	GridGenerator::hyper_shell (triangulation,
				    Point<dim>(), 0.5, 1.0);

	static HyperShellBoundary<dim> boundary;
	triangulation.set_boundary (0, boundary);

	triangulation.refine_global (2);

	break;
      }

      default:
	    Assert (false, ExcNotImplemented());
    }


  const bool do_adaptivity = false;

  if (do_adaptivity)
    {
      setup_dofs(false);

      VectorTools::project (dof_handler,
			    hanging_node_constraints,
			    QGauss<dim>(degree+2),
			    InitialValues<dim>(),
			    old_solution);

      for (unsigned int pre_refinement=0; pre_refinement<4-dim; ++pre_refinement)
	{
	  refine_mesh ();

	  VectorTools::project (dof_handler,
				hanging_node_constraints,
				QGauss<dim>(degree+2),
				InitialValues<dim>(),
				old_solution);
	}
    }
  else
    {
      setup_dofs(true);

      VectorTools::project (dof_handler,
			    hanging_node_constraints,
			    QGauss<dim>(degree+2),
			    InitialValues<dim>(),
			    old_solution);
    }

  timestep_number = 0;
  double time = 0;

  do
    {
      std::cout << "Timestep " << timestep_number
		<< ":  t=" << time
		<< ", dt=" << time_step
                << std::endl;

      std::cout << "   Assembling..." << std::endl;
      assemble_system ();

      std::cout << "   Solving..." << std::endl;
      solve ();

      output_results ();

      time += time_step;
      ++timestep_number;

      old_solution = solution;

      std::cout << std::endl;

      if (do_adaptivity)
	if (timestep_number % 10 == 0)
	  refine_mesh ();
    }
  while (time <= 5);
}



int main ()
{
  try
    {
      deallog.depth_console (0);

      BoussinesqFlowProblem<2> flow_problem(1);
      flow_problem.run ();
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
