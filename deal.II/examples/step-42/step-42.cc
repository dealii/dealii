
/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 1999 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyrightG and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

                                 // @sect3{Include files}

				 // The first few (many?) include
				 // files have already been used in
				 // the previous example, so we will
				 // not explain their meaning here
				 // again.

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vectors.h>
#include <deal.II/numerics/matrices.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <list>
#include <time.h>

				 // This is new, however: in the previous
				 // example we got some unwanted output from
				 // the linear solvers. If we want to suppress
				 // it, we have to include this file and add a
				 // single line somewhere to the program (see
				 // the main() function below for that):
#include <deal.II/base/logstream.h>

				 // The final step, as in previous
				 // programs, is to import all the
				 // deal.II class and function names
				 // into the global namespace:
using namespace dealii;

                                 // @sect3{The <code>Step4</code> class template}

template <int dim> class ConstitutiveLaw;

template <int dim>
class Step4
{
public:
  Step4 ();
  void run ();

private:
  void make_grid ();
  void setup_system();
  void assemble_mass_matrix ();
  void assemble_nl_system (TrilinosWrappers::MPI::Vector &u);
  void residual_nl_system (TrilinosWrappers::MPI::Vector &u,
			   Vector<double>                &sigma_eff_vector);
  void projection_active_set ();
  void dirichlet_constraints ();
  void solve ();
  void solve_newton ();
  void output_results (const std::string& title) const;
  void move_mesh (const TrilinosWrappers::MPI::Vector &_complete_displacement) const;
  void output_results (TrilinosWrappers::MPI::Vector vector, const std::string& title) const;
  void output_results (Vector<double> vector, const std::string& title) const;

  MPI_Comm             mpi_communicator;

  parallel::distributed::Triangulation<dim>   triangulation;

  FESystem<dim>        fe;
  DoFHandler<dim>      dof_handler;

  IndexSet             locally_owned_dofs;
  IndexSet             locally_relevant_dofs;

  int                  n_refinements;
  int                  n_refinements_local;
  unsigned int         number_iterations;
  std::vector<double>  run_time;

  ConstraintMatrix     constraints;
  ConstraintMatrix     constraints_hanging_nodes;
  ConstraintMatrix     constraints_dirichlet_hanging_nodes;

  TrilinosWrappers::SparseMatrix system_matrix_newton;
  TrilinosWrappers::SparseMatrix mass_matrix;

  TrilinosWrappers::MPI::Vector       solution;
  TrilinosWrappers::MPI::Vector       old_solution;
  TrilinosWrappers::MPI::Vector       system_rhs_newton;
  TrilinosWrappers::MPI::Vector       resid_vector;
  TrilinosWrappers::MPI::Vector       diag_mass_matrix_vector;
  IndexSet                            active_set;

  ConditionalOStream pcout;

  TrilinosWrappers::PreconditionAMG::AdditionalData additional_data;
  TrilinosWrappers::PreconditionAMG preconditioner_u;
  TrilinosWrappers::PreconditionAMG preconditioner_t;

  std::auto_ptr<ConstitutiveLaw<dim> > plast_lin_hard;

  double sigma_0;    // Yield stress
  double gamma;      // Parameter for the linear isotropic hardening
  double e_modul;    // E-Modul
  double nu;         // Poisson ratio

  std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionAMG>  Mp_preconditioner;
};

template <int dim>
class ConstitutiveLaw
{
public:
  ConstitutiveLaw (double _E, double _nu, double _sigma_0, double _gamma, MPI_Comm _mpi_communicator, ConditionalOStream _pcout);
  //     ConstitutiveLaw (double mu, double kappa);
  void plast_linear_hardening (SymmetricTensor<4,dim>  &stress_strain_tensor,
			       SymmetricTensor<2,dim>  &strain_tensor,
			       unsigned int 	       &elast_points,
			       unsigned int 	       &plast_points,
			       double                  &sigma_eff,
			       double                  &yield);
  void linearized_plast_linear_hardening (SymmetricTensor<4,dim>  &stress_strain_tensor_linearized,
					  SymmetricTensor<4,dim>  &stress_strain_tensor,
					  SymmetricTensor<2,dim>  &strain_tensor);
  inline SymmetricTensor<2,dim> get_strain (const FEValues<dim> &fe_values,
					    const unsigned int  shape_func,
					    const unsigned int  q_point) const;

private:
  SymmetricTensor<4,dim>  stress_strain_tensor_mu;
  SymmetricTensor<4,dim>  stress_strain_tensor_kappa;
  double E;
  double nu;
  double sigma_0;
  double gamma;
  double mu;
  double kappa;
  MPI_Comm mpi_communicator;
  ConditionalOStream pcout;
};

template <int dim>
ConstitutiveLaw<dim>::ConstitutiveLaw(double _E, double _nu, double _sigma_0, double _gamma, MPI_Comm _mpi_communicator, ConditionalOStream _pcout)
 :E (_E),
  nu (_nu),
  sigma_0 (_sigma_0),
  gamma (_gamma),
  mpi_communicator (_mpi_communicator),
  pcout (_pcout)
{
  mu = E/(2*(1+nu));
  kappa = E/(3*(1-2*nu));
  pcout<< "-----> mu = " << mu << ", kappa = " << kappa <<std::endl;
  stress_strain_tensor_kappa = kappa*outer_product(unit_symmetric_tensor<dim>(), unit_symmetric_tensor<dim>());
  stress_strain_tensor_mu = 2*mu*(identity_tensor<dim>() - outer_product(unit_symmetric_tensor<dim>(), unit_symmetric_tensor<dim>())/3.0);
}

template <int dim>
inline
SymmetricTensor<2,dim> ConstitutiveLaw<dim>::get_strain (const FEValues<dim> &fe_values,
						     const unsigned int   shape_func,
						     const unsigned int   q_point) const
{
  const FEValuesExtractors::Vector displacement (0);
  SymmetricTensor<2,dim> tmp;

  tmp = fe_values[displacement].symmetric_gradient (shape_func,q_point);

  return tmp;
}

template <int dim>
void ConstitutiveLaw<dim>::plast_linear_hardening (SymmetricTensor<4,dim>  &stress_strain_tensor,
						   SymmetricTensor<2,dim>  &strain_tensor,
						   unsigned int            &elast_points,
						   unsigned int            &plast_points,
						   double                  &sigma_eff,
						   double                  &yield)
{
  // Plane strain
  if (dim == 3)
  {
    SymmetricTensor<2,dim> stress_tensor;
    stress_tensor = (stress_strain_tensor_kappa + stress_strain_tensor_mu)*strain_tensor;
    double tmp = E/((1+nu)*(1-2*nu));
    double stress_tensor_33 = 0.0;//tmp*(strain_tensor[0][0] + strain_tensor[1][1])*nu;

    SymmetricTensor<2,dim> deviator_stress_tensor = deviator(stress_tensor);

    double deviator_stress_tensor_norm = deviator_stress_tensor.norm ();
    deviator_stress_tensor_norm = std::sqrt (deviator_stress_tensor_norm*deviator_stress_tensor_norm +
				  stress_tensor_33*stress_tensor_33);

    yield = 0;
    stress_strain_tensor = stress_strain_tensor_mu;
    double beta = 1.0;
    if (deviator_stress_tensor_norm >= sigma_0)
    {
      beta = (sigma_0 + gamma)/deviator_stress_tensor_norm;
      stress_strain_tensor *= beta;
      yield = 1;
      plast_points += 1;
    }
    else
      elast_points += 1;

//     std::cout<< beta <<std::endl;
    stress_strain_tensor += stress_strain_tensor_kappa;

    sigma_eff = beta * deviator_stress_tensor_norm;
  }
}

template <int dim>
void ConstitutiveLaw<dim>::linearized_plast_linear_hardening (SymmetricTensor<4,dim>  &stress_strain_tensor_linearized,
							      SymmetricTensor<4,dim>  &stress_strain_tensor,
							      SymmetricTensor<2,dim>  &strain_tensor)
{
  // Plane strains
  if (dim == 3)
  {
    SymmetricTensor<2,dim> stress_tensor;
    stress_tensor = (stress_strain_tensor_kappa + stress_strain_tensor_mu)*strain_tensor;
    double tmp = E/((1+nu)*(1-2*nu));
    double stress_tensor_33 = 0.0;//tmp*(strain_tensor[0][0] + strain_tensor[1][1])*nu;

    stress_strain_tensor = stress_strain_tensor_mu;
    stress_strain_tensor_linearized = stress_strain_tensor_mu;

    SymmetricTensor<2,dim> deviator_stress_tensor = deviator(stress_tensor);

    double deviator_stress_tensor_norm = deviator_stress_tensor.norm ();
    deviator_stress_tensor_norm = std::sqrt (deviator_stress_tensor_norm*deviator_stress_tensor_norm + stress_tensor_33*stress_tensor_33);
    double beta = 1.0;
    if (deviator_stress_tensor_norm >= sigma_0)
    {
      beta = (sigma_0 + gamma)/deviator_stress_tensor_norm;
      stress_strain_tensor *= beta;
      stress_strain_tensor_linearized *= beta;
      deviator_stress_tensor /= deviator_stress_tensor_norm;
      stress_strain_tensor_linearized -= beta*2*mu*outer_product(deviator_stress_tensor, deviator_stress_tensor);
    }

    stress_strain_tensor += stress_strain_tensor_kappa;
    stress_strain_tensor_linearized += stress_strain_tensor_kappa;
  }
}

namespace EquationData
{
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide () : Function<dim>(dim) {}

    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
			       Vector<double>   &values) const;
  };

  template <int dim>
  double RightHandSide<dim>::value (const Point<dim> &p,
				    const unsigned int component) const
  {
    double return_value = 0.0;

    if (component == 0)
      return_value = 0.0;
    if (component == 1)
      return_value = 0.0;
    if (component == 2)
      // if ((p(0)-0.5)*(p(0)-0.5)+(p(1)-0.5)*(p(1)-0.5) < 0.2)
      // 	return_value = -5000;
      // else
      return_value = 0.0;
    // for (unsigned int i=0; i<dim; ++i)
    //   return_value += 4*std::pow(p(i), 4);

    return return_value;
  }

  template <int dim>
  void RightHandSide<dim>::vector_value (const Point<dim> &p,
					 Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = RightHandSide<dim>::value (p, c);
  }


  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues () : Function<dim>(dim) {};

    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
			       Vector<double>   &values) const;
  };

  template <int dim>
  double BoundaryValues<dim>::value (const Point<dim> &p,
				     const unsigned int component) const
  {
    double return_value = 0;

    if (component == 0)
      return_value = 0.0;
    if (component == 1)
      return_value = 0.0;
    if (component == 2)
      return_value = 0.0;

    return return_value;
  }

  template <int dim>
  void BoundaryValues<dim>::vector_value (const Point<dim> &p,
					   Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = BoundaryValues<dim>::value (p, c);
  }


  template <int dim>
  class Obstacle : public Function<dim>
  {
  public:
    Obstacle () : Function<dim>(dim) {};

    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
			       Vector<double>   &values) const;
  };

  template <int dim>
  double Obstacle<dim>::value (const Point<dim> &p,
			       const unsigned int component) const
  {
    double R = 0.03;
    double return_value = 0.0;
    if (component == 0)
      return_value = p(0);
    if (component == 1)
      return_value = p(1);
    if (component == 2)
      {
	// double hz = 0.98;
	// double position_x = 0.5;
	// double alpha = 12.0;
	// double s_x = 0.5039649116;
	// double s_y = hz + 0.00026316298;
	// if (p(0) > position_x - R && p(0) < s_x)
	//   {
	//     return_value = -sqrt(R*R - (p(0)-position_x)*(p(0)-position_x)) + hz + R;
	//   }
	// else if (p(0) >= s_x)
	//   {
	//     return_value = 12.0/90.0*p(0) + (s_y - alpha/90.0*s_x);
	//   }
	// else
	//   return_value = 1e+10;

	// Hindernis Dortmund
	// double x1 = p(0);
	// double x2 = p(1);
	// if (((x2-0.5)*(x2-0.5)+(x1-0.5)*(x1-0.5)<=0.3*0.3)&&((x2-0.5)*(x2-0.5)+(x1-1.0)*(x1-1.0)>=0.4*0.4)&&((x2-0.5)*(x2-0.5)+x1*x1>=0.4*0.4))
	//   return_value = 0.999;
	// else
	//   return_value = 1e+10;

	// Hindernis Werkzeug TKSE
	// double shift_walze_x = 0.0;
	// double shift_walze_y = 0.0;
	// return_value = 0.032 + data->dicke - input_copy->mikro_height (p(0) + shift_walze_x, p(1) + shift_walze_y, p(2));

	// Ball with radius R
	double R = 0.5;
	if (std::pow ((p(0)-1.0/2.0), 2) + std::pow ((p(1)-1.0/2.0), 2) < R*R)
	  return_value = 1.0 + R - 0.001 - sqrt (R*R  - std::pow ((p(0)-1.0/2.0), 2)
						 - std::pow ((p(1)-1.0/2.0), 2));
	else
	  return_value = 1e+5;
      }
    return return_value;

    // return 1e+10;//0.98;
  }

  template <int dim>
  void Obstacle<dim>::vector_value (const Point<dim> &p,
			            Vector<double> &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = Obstacle<dim>::value (p, c);
  }
}


                                 // @sect3{Implementation of the <code>Step4</code> class}

                                 // Next for the implementation of the class
                                 // template that makes use of the functions
                                 // above. As before, we will write everything

template <int dim>
Step4<dim>::Step4 ()
  :
  mpi_communicator (MPI_COMM_WORLD),
  triangulation (mpi_communicator),
  fe (FE_Q<dim>(1), dim),
  dof_handler (triangulation),
  pcout (std::cout,
	 (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
  sigma_0 (400),
  gamma (1.e-2),
  e_modul (2.e5),
  nu (0.3)
{
  // double _E, double _nu, double _sigma_0, double _gamma
  plast_lin_hard.reset (new ConstitutiveLaw<dim> (e_modul, nu, sigma_0, gamma, mpi_communicator, pcout));
}

template <int dim>
void Step4<dim>::make_grid ()
{
  std::vector<unsigned int> repet(3);
  repet[0] = 1;//20;
  repet[1] = 1;
  repet[2] = 1;

  Point<dim> p1 (0,0,0);
  Point<dim> p2 (1.0, 1.0, 1.0);
  GridGenerator::subdivided_hyper_rectangle (triangulation, repet, p1, p2);

  Triangulation<3>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();

  /* boundary_indicators:
            _______
           /  9    /|
          /______ / |
        8|       | 8|
         |   8   | /
         |_______|/
             6
   */

  for (; cell!=endc; ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
      {
	if (cell->face (face)->center ()[2] == p2(2))
	  cell->face (face)->set_boundary_indicator (9);
	if (cell->face (face)->center ()[0] == p1(0) ||
	    cell->face (face)->center ()[0] == p2(0) ||
	    cell->face (face)->center ()[1] == p1(1) ||
	    cell->face (face)->center ()[1] == p2(1))
	  cell->face (face)->set_boundary_indicator (8);
	if (cell->face (face)->center ()[2] == p1(2))
	  cell->face (face)->set_boundary_indicator (6);
      }

  n_refinements = 3;
  n_refinements_local = 3;
  triangulation.refine_global (n_refinements);

  // Lokale Verfeinerung des Gitters
  for (int step=0; step<n_refinements_local; ++step)
    {
      cell = triangulation.begin_active();  // Iterator ueber alle Zellen

      double hlp_refinement = 0;
      hlp_refinement = pow((double)(step)/(n_refinements_local),4.0);
      pcout<< "Verfeinerungsfaktor: " << hlp_refinement <<std::endl;

      for (; cell!=endc; ++cell)
         for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
         {
	   if (cell->face (face)->at_boundary()
	       && cell->face (face)->boundary_indicator () == 9)
	     {
	       cell->set_refine_flag ();
	       break;
	     }
	   else if (cell->level () == n_refinements + n_refinements_local - 1)
	     {
	       cell->set_refine_flag ();
	       break;
	     }
        };
      triangulation.execute_coarsening_and_refinement ();
    };
}

template <int dim>
void Step4<dim>::setup_system ()
{
  // setup dofs
  {
    dof_handler.distribute_dofs (fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    locally_relevant_dofs.clear();
    DoFTools::extract_locally_relevant_dofs (dof_handler,
					     locally_relevant_dofs);
  }

  // setup hanging nodes and dirichlet constraints
  {
    // constraints_hanging_nodes.clear ();
    constraints_hanging_nodes.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler,
					     constraints_hanging_nodes);
    constraints_hanging_nodes.close ();

    pcout << "Number of active cells: "
	       << triangulation.n_active_cells()
	       << std::endl
	       << "Total number of cells: "
	       << triangulation.n_cells()
	       << std::endl
	       << "Number of degrees of freedom: "
	       << dof_handler.n_dofs ()
	       << std::endl;

    dirichlet_constraints ();
  }

  // Initialzation for matrices and vectors
  {
    solution.reinit (locally_relevant_dofs, mpi_communicator);
    system_rhs_newton.reinit (locally_owned_dofs, mpi_communicator);
    old_solution.reinit (system_rhs_newton);
    resid_vector.reinit (system_rhs_newton);
    diag_mass_matrix_vector.reinit (system_rhs_newton);
    active_set.set_size (locally_relevant_dofs.size ());
  }

  // setup sparsity pattern
  {
    TrilinosWrappers::SparsityPattern sp (locally_owned_dofs,
					  mpi_communicator);

    DoFTools::make_sparsity_pattern (dof_handler, sp, constraints_dirichlet_hanging_nodes, false,
				     Utilities::MPI::this_mpi_process(mpi_communicator));

    sp.compress();

    system_matrix_newton.reinit (sp);

    mass_matrix.reinit (sp);
  }

  assemble_mass_matrix ();
  const unsigned int
    start = (system_rhs_newton.local_range().first),
    end   = (system_rhs_newton.local_range().second);
  for (unsigned int j=0; j<end; j++)
    diag_mass_matrix_vector (j) = mass_matrix.diag_element (j);
  number_iterations = 0;
}

template <int dim>
void Step4<dim>::assemble_mass_matrix ()
{
  QTrapez<dim-1>  face_quadrature_formula;

  FEFaceValues<dim> fe_values_face (fe, face_quadrature_formula,
			              update_values   | update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell      = fe.dofs_per_cell;
  const unsigned int   dofs_per_face      = fe.dofs_per_face;
  const unsigned int   n_face_q_points    = face_quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  const FEValuesExtractors::Vector displacement (0);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  for (; cell!=endc; ++cell)
    if (cell->is_locally_owned())
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	if (cell->face (face)->at_boundary()
	    && cell->face (face)->boundary_indicator () == 9)
	  {
	    fe_values_face.reinit (cell, face);
	    cell_matrix = 0;

	    for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		cell_matrix(i,i) += (fe_values_face[displacement].value (i, q_point) *
				     fe_values_face[displacement].value (i, q_point) *
				     fe_values_face.JxW (q_point));

	    cell->get_dof_indices (local_dof_indices);

	    constraints_dirichlet_hanging_nodes.distribute_local_to_global (cell_matrix,
									    local_dof_indices,
									    mass_matrix);
	  }

  mass_matrix.compress ();
}

template <int dim>
void Step4<dim>::assemble_nl_system (TrilinosWrappers::MPI::Vector &u)
{
  QGauss<dim>  quadrature_formula(2);
  QGauss<dim-1>  face_quadrature_formula(2);

  FEValues<dim> fe_values (fe, quadrature_formula,
			   UpdateFlags(update_values    |
				       update_gradients |
				       update_q_points  |
				       update_JxW_values));

  FEFaceValues<dim> fe_values_face (fe, face_quadrature_formula,
				    update_values   | update_quadrature_points |
				    update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size ();
  const unsigned int   n_face_q_points = face_quadrature_formula.size();

  const EquationData::RightHandSide<dim> right_hand_side;
  std::vector<Vector<double> > right_hand_side_values (n_q_points,
                                                       Vector<double>(dim));
  std::vector<Vector<double> > right_hand_side_values_face (n_face_q_points,
							    Vector<double>(dim));

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
						 endc = dof_handler.end();

  const FEValuesExtractors::Vector displacement (0);

  TrilinosWrappers::MPI::Vector   test_rhs(solution);
  const double kappa = 1.0;
  for (; cell!=endc; ++cell)
    if (cell->is_locally_owned())
      {
	fe_values.reinit (cell);
	cell_matrix = 0;
	cell_rhs = 0;

	right_hand_side.vector_value_list (fe_values.get_quadrature_points(),
					   right_hand_side_values);

	std::vector<SymmetricTensor<2,dim> > strain_tensor (n_q_points);
	fe_values[displacement].get_function_symmetric_gradients (u, strain_tensor);

	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  {
	    SymmetricTensor<4,dim> stress_strain_tensor_linearized;
	    SymmetricTensor<4,dim> stress_strain_tensor;
	    SymmetricTensor<2,dim> stress_tensor;

	    plast_lin_hard->linearized_plast_linear_hardening (stress_strain_tensor_linearized,
							       stress_strain_tensor,
							       strain_tensor[q_point]);

	    //   	if (q_point == 0)
	    //  	std::cout<< stress_strain_tensor_linearized <<std::endl;
	    //  	std::cout<< stress_strain_tensor <<std::endl;
	    for (unsigned int i=0; i<dofs_per_cell; ++i)
	      {
		stress_tensor = stress_strain_tensor_linearized * plast_lin_hard->get_strain(fe_values, i, q_point);

		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  {
		    cell_matrix(i,j) += (stress_tensor *
					 plast_lin_hard->get_strain(fe_values, j, q_point) *
					 fe_values.JxW (q_point));
		  }

		// the linearized part a(v^i;v^i,v) of the rhs
		cell_rhs(i) += (stress_tensor *
				strain_tensor[q_point] *
				fe_values.JxW (q_point));

		// the residual part a(v^i;v) of the rhs
		cell_rhs(i) -= (strain_tensor[q_point] * stress_strain_tensor *
				plast_lin_hard->get_strain(fe_values, i, q_point) *
				fe_values.JxW (q_point));

		// the residual part F(v) of the rhs
		Tensor<1,dim> rhs_values;
		rhs_values = 0;
		cell_rhs(i) += (fe_values[displacement].value (i, q_point) *
				rhs_values *
				fe_values.JxW (q_point));
	      }
	  }

	for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	  {
	    if (cell->face (face)->at_boundary()
		&& cell->face (face)->boundary_indicator () == 9)
	      {
		fe_values_face.reinit (cell, face);

		right_hand_side.vector_value_list (fe_values_face.get_quadrature_points(),
						   right_hand_side_values_face);

		for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
		  {
		    Tensor<1,dim> rhs_values;
		    rhs_values = 0;
		    for (unsigned int i=0; i<dofs_per_cell; ++i)
		      cell_rhs(i) += (fe_values_face[displacement].value (i, q_point) *
				      rhs_values *
				      fe_values_face.JxW (q_point));
		  }
	      }
	  }

	cell->get_dof_indices (local_dof_indices);
	constraints.distribute_local_to_global (cell_matrix, cell_rhs,
						local_dof_indices,
						system_matrix_newton, system_rhs_newton, true);
      };

  system_matrix_newton.compress ();
  system_rhs_newton.compress ();
}

template <int dim>
void Step4<dim>::residual_nl_system (TrilinosWrappers::MPI::Vector &u,
				     Vector<double>                &sigma_eff_vector)
{
  QGauss<dim>  quadrature_formula(2);
  QGauss<dim-1> face_quadrature_formula(2);

  FEValues<dim> fe_values (fe, quadrature_formula,
			   UpdateFlags(update_values    |
				       update_gradients |
				       update_q_points  |
				       update_JxW_values));

  FEFaceValues<dim> fe_values_face (fe, face_quadrature_formula,
				    update_values   | update_quadrature_points |
				    update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size ();
  const unsigned int   n_face_q_points = face_quadrature_formula.size();

  const EquationData::RightHandSide<dim> right_hand_side;
  std::vector<Vector<double> > right_hand_side_values (n_q_points,
                                                       Vector<double>(dim));
  std::vector<Vector<double> > right_hand_side_values_face (n_face_q_points,
							    Vector<double>(dim));

  Vector<double>       cell_rhs (dofs_per_cell);
  Vector<double>       cell_sigma_eff (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  const FEValuesExtractors::Vector displacement (0);

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
						 endc = dof_handler.end();

  unsigned int elast_points = 0;
  unsigned int plast_points = 0;
  double       sigma_eff = 0;
  double       yield = 0;
  unsigned int cell_number = 0;
  for (; cell!=endc; ++cell)
    if (cell->is_locally_owned())
      {
	fe_values.reinit (cell);
	cell_rhs = 0;

	right_hand_side.vector_value_list (fe_values.get_quadrature_points(),
					   right_hand_side_values);

	std::vector<SymmetricTensor<2,dim> > strain_tensor (n_q_points);
	fe_values[displacement].get_function_symmetric_gradients (u, strain_tensor);

	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  {
	    SymmetricTensor<4,dim> stress_strain_tensor;
	    SymmetricTensor<2,dim> stress_tensor;

	    plast_lin_hard->plast_linear_hardening (stress_strain_tensor, strain_tensor[q_point],
						    elast_points, plast_points, sigma_eff, yield);

	    // sigma_eff_vector (cell_number) += sigma_eff;
	    sigma_eff_vector (cell_number) += yield;

	    /*	if (q_point == 0)
		std::cout<< stress_strain_tensor <<std::endl;*/
	    for (unsigned int i=0; i<dofs_per_cell; ++i)
	      {
		cell_rhs(i) -= (strain_tensor[q_point] * stress_strain_tensor * //(stress_tensor) *
				plast_lin_hard->get_strain(fe_values, i, q_point) *
				fe_values.JxW (q_point));

		Tensor<1,dim> rhs_values;
		rhs_values = 0;
		cell_rhs(i) += ((fe_values[displacement].value (i, q_point) *
				 rhs_values) *
				fe_values.JxW (q_point));
	      };
	  };

	for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	  {
	    if (cell->face (face)->at_boundary()
		&& cell->face (face)->boundary_indicator () == 9)
	      {
		fe_values_face.reinit (cell, face);

		right_hand_side.vector_value_list (fe_values_face.get_quadrature_points(),
						   right_hand_side_values_face);

		for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
		  {
		    Tensor<1,dim> rhs_values;
		    rhs_values = 0;
		    for (unsigned int i=0; i<dofs_per_cell; ++i)
		      cell_rhs(i) += (fe_values_face[displacement].value (i, q_point) *
				      rhs_values *
				      fe_values_face.JxW (q_point));
		  }
	      }
	  }

	cell->get_dof_indices (local_dof_indices);
	constraints_dirichlet_hanging_nodes.distribute_local_to_global (cell_rhs,
									local_dof_indices,
									system_rhs_newton);

	sigma_eff_vector(cell_number) /= n_q_points;
	cell_number += 1;
      };

  system_rhs_newton.compress ();

  unsigned int sum_elast_points = Utilities::MPI::sum(elast_points, mpi_communicator);
  unsigned int sum_plast_points = Utilities::MPI::sum(plast_points, mpi_communicator);
  pcout<< "Elast-Points = " << sum_elast_points <<std::endl;
  pcout<< "Plast-Points = " << sum_plast_points <<std::endl;
}

                                 // @sect4{Step4::projection_active_set}

				 // Projection and updating of the active set
                                 // for the dofs which penetrates the obstacle.
template <int dim>
void Step4<dim>::projection_active_set ()
{
  const EquationData::Obstacle<dim>     obstacle;
  std::vector<bool>                     vertex_touched (dof_handler.n_dofs (), false);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();

  TrilinosWrappers::MPI::Vector         distributed_solution (system_rhs_newton);
  distributed_solution = solution;
  TrilinosWrappers::MPI::Vector         lambda (solution);
  lambda = resid_vector;
  TrilinosWrappers::MPI::Vector         diag_mass_matrix_vector_relevant (solution);
  diag_mass_matrix_vector_relevant = diag_mass_matrix_vector;

  constraints.reinit(locally_relevant_dofs);
  active_set.clear ();
  const double c = 100.0*e_modul;

  for (; cell!=endc; ++cell)
    if (cell->is_locally_owned())
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	if (cell->face (face)->boundary_indicator () == 9)
	  for (unsigned int v=0; v<GeometryInfo<dim-1>::vertices_per_cell; ++v)
	    {
	      unsigned int index_z = cell->face (face)->vertex_dof_index (v,2);

	      if (vertex_touched[cell->face (face)->vertex_index(v)] == false)
		vertex_touched[cell->face (face)->vertex_index(v)] = true;
	      else
		continue;

	      // the local row where
	      Point<dim> point (cell->face (face)->vertex (v)[0],/* + solution (index_x),*/
				cell->face (face)->vertex (v)[1],
				cell->face (face)->vertex (v)[2]);

	      double obstacle_value = obstacle.value (point, 2);
	      double solution_index_z = solution (index_z);
	      double gap = obstacle_value - point (2);

	      if (lambda (index_z) +
		  c*diag_mass_matrix_vector_relevant (index_z)*(solution_index_z - gap) > 0)
		{
		  constraints.add_line (index_z);
		  constraints.set_inhomogeneity (index_z, gap);

		  distributed_solution (index_z) = gap;

		  if (locally_relevant_dofs.is_element (index_z))
		    active_set.add_index (index_z);

		  // std::cout<< index_z << ", "
		  // 	   << "Error: " << lambda (index_z) +
		  //   diag_mass_matrix_vector_relevant (index_z)*c*(solution_index_z - gap)
		  // 	   << ", " << lambda (index_z)
		  // 	   << ", " << diag_mass_matrix_vector_relevant (index_z)
		  // 	   << ", " << obstacle_value
		  // 	   << ", " << solution_index_z
		  // 	   <<std::endl;
		}
	    }

  distributed_solution.compress(Insert);

//TODO: since the active sets on different processors also store "ghost" elements,
// i.e. constraints for DoFs that are locally_relevant but not locally_owned,
// the following sum doesn't make much sense. We need to count the
// *unique* elements
  unsigned int sum_contact_constraints = Utilities::MPI::sum(active_set.n_elements (), mpi_communicator);
  pcout << "Number of Contact-Constaints: " << sum_contact_constraints <<std::endl;

  solution = distributed_solution;

  constraints.close ();

  const ConstraintMatrix::MergeConflictBehavior
    merge_conflict_behavior = ConstraintMatrix::left_object_wins;
  constraints.merge (constraints_dirichlet_hanging_nodes, merge_conflict_behavior);
}

template <int dim>
void Step4<dim>::dirichlet_constraints ()
{
  /* boundary_indicators:
            _______
           /  9    /|
          /______ / |
        8|       | 8|
         |   8   | /
         |_______|/
             6
   */

  constraints_dirichlet_hanging_nodes.reinit (locally_relevant_dofs);
  constraints_dirichlet_hanging_nodes.merge (constraints_hanging_nodes);

  std::vector<bool> component_mask (dim, true);
  component_mask[0] = true;
  component_mask[1] = true;
  component_mask[2] = true;
  VectorTools::interpolate_boundary_values (dof_handler,
					    6,
					    EquationData::BoundaryValues<dim>(),
					    constraints_dirichlet_hanging_nodes,
					    component_mask);

  component_mask[0] = true;
  component_mask[1] = true;
  component_mask[2] = false;
  VectorTools::interpolate_boundary_values (dof_handler,
  					    8,
  					    EquationData::BoundaryValues<dim>(),
  					    constraints_dirichlet_hanging_nodes,
  					    component_mask);
  constraints_dirichlet_hanging_nodes.close ();
}

template <int dim>
void Step4<dim>::solve ()
{
  ReductionControl                 reduction_control (10000, 1e-15, 1e-4);

  TrilinosWrappers::MPI::Vector    distributed_solution (system_rhs_newton);
  distributed_solution = solution;

  constraints_hanging_nodes.set_zero (distributed_solution);

  // Solving iterative
  SolverCG<TrilinosWrappers::MPI::Vector>
    solver (reduction_control, mpi_communicator);

  preconditioner_u.initialize (system_matrix_newton, additional_data);

  solver.solve (system_matrix_newton, distributed_solution, system_rhs_newton, preconditioner_u);
  pcout << "Initial error: " << reduction_control.initial_value() <<std::endl;
  pcout << "   " << reduction_control.last_step()
  	    << " CG iterations needed to obtain convergence with an error: "
  	    <<  reduction_control.last_value()
  	    << std::endl;

  number_iterations += reduction_control.last_step();

  constraints.distribute (distributed_solution);

  solution = distributed_solution;
}

template <int dim>
void Step4<dim>::solve_newton ()
{
  double                         resid=0;
  double                         resid_old=100000;
  TrilinosWrappers::MPI::Vector  res (system_rhs_newton);
  TrilinosWrappers::MPI::Vector  tmp_vector (system_rhs_newton);
  clock_t                        start, end;

  std::vector<std::vector<bool> > constant_modes;
  std::vector<bool>  components (dim,true);
  components[dim] = false;
  DoFTools::extract_constant_modes (dof_handler, components,
				    constant_modes);

  additional_data.elliptic = true;
  additional_data.n_cycles = 1;
  additional_data.w_cycle = false;
  additional_data.output_details = false;
  additional_data.smoother_sweeps = 2;
  additional_data.aggregation_threshold = 1e-2;

  IndexSet                            active_set_old (active_set);
  Vector<double>                      sigma_eff_vector;
  sigma_eff_vector.reinit (triangulation.n_active_cells());
  unsigned int j = 0;
  unsigned int number_assemble_system = 0;
  for (; j<=100;j++)
    {
      pcout<< " " <<std::endl;
      pcout<< j << ". Iteration of the inexact Newton-method." <<std::endl;
      pcout<< "Update of active set" <<std::endl;
      projection_active_set ();

      pcout<< "Assembling ... " <<std::endl;
      start = clock();
      system_matrix_newton = 0;
      system_rhs_newton = 0;
      assemble_nl_system (solution);  //compute Newton-Matrix
      end = clock();
      run_time[1] += (double)(end-start)/CLOCKS_PER_SEC;

      number_assemble_system += 1;

      start = clock();
      solve ();
      end = clock();
      run_time[2] += (double)(end-start)/CLOCKS_PER_SEC;

      TrilinosWrappers::MPI::Vector    distributed_solution (system_rhs_newton);
      distributed_solution = solution;

      int damped = 0;
      tmp_vector = old_solution;
      double a = 0;
      for (unsigned int i=0; (i<10)&&(!damped); i++)
	{
	  a=pow(0.5,i);
	  old_solution = tmp_vector;
	  old_solution.sadd(1-a,a, distributed_solution);

	  start = clock();
	  system_rhs_newton = 0;
	  sigma_eff_vector = 0;
	  solution = old_solution;
	  residual_nl_system (solution, sigma_eff_vector);
	  res = system_rhs_newton;

	  const unsigned int
	    start_res     = (res.local_range().first),
	    end_res       = (res.local_range().second);
	  for (unsigned int n=start_res; n<end_res; ++n)
	    if (constraints.is_inhomogeneously_constrained (n))
	      {
		// pcout<< i << ". " << constraints.get_inhomogeneity (n)
		// 	 << ". " << res (n)
		// 	 << ", start = " << start_res
		// 	 << ", end = " << end_res
		// 	 <<std::endl;
		res(n) = 0;
	      }

	  resid = res.l2_norm ();
	  pcout<< "Residual: " << resid <<std::endl;

	  if (resid<resid_old)
	    {
	      pcout<< "Newton-damping parameter alpha = " << a <<std::endl;
	      damped=1;
	    }
	  end = clock();
	  run_time[3] = (double)(end-start)/CLOCKS_PER_SEC;
	}

      if (resid<1e-8)
	{
	  pcout<< "Inexact Newton-method stopped with residual = " << resid <<std::endl;
	  pcout<< "Number of Assembling systems = " << number_assemble_system <<std::endl;
	  break;
	}
      resid_old=resid;

      resid_vector = system_rhs_newton;

      if (active_set == active_set_old && resid < 1e-10)
	break;
      active_set_old = active_set;
    } // End of active-set-loop

  start = clock();
  pcout<< "Creating output." <<std::endl;
  std::ostringstream filename_solution;
  filename_solution << "solution";
  // filename_solution << "solution_";
  // filename_solution << k;
  output_results (filename_solution.str ());
  // output_results (sigma_eff_vector, "sigma_eff");
  end = clock();
  run_time[4] = (double)(end-start)/CLOCKS_PER_SEC;

  pcout<< "Number of Solver-Iterations = " << number_iterations <<std::endl;

  pcout<< "%%%%%% Rechenzeit make grid and setup = " << run_time[0] <<std::endl;
  pcout<< "%%%%%% Rechenzeit assemble system = " << run_time[1] <<std::endl;
  pcout<< "%%%%%% Rechenzeit solve system = " << run_time[2] <<std::endl;
  pcout<< "%%%%%% Rechenzeit error and lambda = " << run_time[3] <<std::endl;
  pcout<< "%%%%%% Rechenzeit output = " << run_time[4] <<std::endl;
}

template <int dim>
void Step4<dim>::output_results (const std::string& title) const
{
  move_mesh (solution);

  TrilinosWrappers::MPI::Vector         lambda (solution);
  lambda = resid_vector;

  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);

  const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
  data_out.add_data_vector (solution, std::vector<std::string>(dim, "Displacement"),
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);
  data_out.add_data_vector (lambda, std::vector<std::string>(dim, "Residual"),
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);
  data_out.add_data_vector (active_set, std::vector<std::string>(dim, "ActiveSet"),
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);

  Vector<float> subdomain (triangulation.n_active_cells());
  for (unsigned int i=0; i<subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector (subdomain, "subdomain");

  data_out.build_patches ();

  const std::string filename = (title + "-" +
  				Utilities::int_to_string
  				(triangulation.locally_owned_subdomain(), 4));

  std::ofstream output_vtu ((filename + ".vtu").c_str ());
  data_out.write_vtu (output_vtu);

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::vector<std::string> filenames;
      for (unsigned int i=0;
  	   i<Utilities::MPI::n_mpi_processes(mpi_communicator);
  	   ++i)
  	filenames.push_back ("solution-" +
  			     Utilities::int_to_string (i, 4) +
  			     ".vtu");

      std::ofstream master_output ((filename + ".pvtu").c_str());
      data_out.write_pvtu_record (master_output, filenames);
    }

  TrilinosWrappers::MPI::Vector  tmp (solution);
  tmp *= -1;
  move_mesh (tmp);
}

template <int dim>
void Step4<dim>::move_mesh (const TrilinosWrappers::MPI::Vector &_complete_displacement) const
{
  pcout<< "Moving mesh." <<std::endl;

  std::vector<bool> vertex_touched (triangulation.n_vertices(),
				    false);

  for (typename DoFHandler<dim>::active_cell_iterator
	 cell = dof_handler.begin_active ();
       cell != dof_handler.end(); ++cell)
    if (cell->is_locally_owned())
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
	{
	  if (vertex_touched[cell->vertex_index(v)] == false)
	    {
	      vertex_touched[cell->vertex_index(v)] = true;

	      Point<dim> vertex_displacement;
	      for (unsigned int d=0; d<dim; ++d)
		{
		  if (_complete_displacement(cell->vertex_dof_index(v,d)) != 0)
		    vertex_displacement[d]
		      = _complete_displacement(cell->vertex_dof_index(v,d));
		}

	      cell->vertex(v) += vertex_displacement;
	    }
	}
}

template <int dim>
void Step4<dim>::output_results (TrilinosWrappers::MPI::Vector vector,
				 const std::string& title) const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (vector, "vector_to_plot");

  data_out.build_patches ();

  std::ofstream output_vtk (dim == 2 ?
			    (title + ".vtk").c_str () :
			    (title + ".vtk").c_str ());
  data_out.write_vtk (output_vtk);
}

template <int dim>
void Step4<dim>::output_results (Vector<double> vector, const std::string& title) const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (vector, "vector_to_plot");

  data_out.build_patches ();

  std::ofstream output_vtk (dim == 2 ?
			    (title + ".vtk").c_str () :
			    (title + ".vtk").c_str ());
  data_out.write_vtk (output_vtk);
}

template <int dim>
void Step4<dim>::run ()
{
  pcout << "Solving problem in " << dim << " space dimensions." << std::endl;

  run_time.resize (5);

  clock_t     start, end;

  start = clock();
  make_grid();
  //  mesh_surface ();
  setup_system ();
  end = clock();
  run_time[0] = (double)(end-start)/CLOCKS_PER_SEC;

  solve_newton ();
}


                                 // @sect3{The <code>main</code> function}

				 // And this is the main function. It also
				 // looks mostly like in step-3, but if you
				 // look at the code below, note how we first
				 // create a variable of type
				 // <code>Step4@<2@></code> (forcing
				 // the compiler to compile the class template
				 // with <code>dim</code> replaced by
				 // <code>2</code>) and run a 2d simulation,
				 // and then we do the whole thing over in 3d.
				 //
				 // In practice, this is probably not what you
				 // would do very frequently (you probably
				 // either want to solve a 2d problem, or one
				 // in 3d, but not both at the same
				 // time). However, it demonstrates the
				 // mechanism by which we can simply change
				 // which dimension we want in a single place,
				 // and thereby force the compiler to
				 // recompile the dimension independent class
				 // templates for the dimension we
				 // request. The emphasis here lies on the
				 // fact that we only need to change a single
				 // place. This makes it rather trivial to
				 // debug the program in 2d where computations
				 // are fast, and then switch a single place
				 // to a 3 to run the much more computing
				 // intensive program in 3d for `real'
				 // computations.
				 //
				 // Each of the two blocks is enclosed in
				 // braces to make sure that the
				 // <code>laplace_problem_2d</code> variable
				 // goes out of scope (and releases the memory
				 // it holds) before we move on to allocate
				 // memory for the 3d case. Without the
				 // additional braces, the
				 // <code>laplace_problem_2d</code> variable
				 // would only be destroyed at the end of the
				 // function, i.e. after running the 3d
				 // problem, and would needlessly hog memory
				 // while the 3d run could actually use it.
                                 //
                                 // Finally, the first line of the function is
                                 // used to suppress some output.  Remember
                                 // that in the previous example, we had the
                                 // output from the linear solvers about the
                                 // starting residual and the number of the
                                 // iteration where convergence was
                                 // detected. This can be suppressed through
                                 // the <code>deallog.depth_console(0)</code>
                                 // call.
                                 //
                                 // The rationale here is the following: the
                                 // deallog (i.e. deal-log, not de-allog)
                                 // variable represents a stream to which some
                                 // parts of the library write output. It
                                 // redirects this output to the console and
                                 // if required to a file. The output is
                                 // nested in a way so that each function can
                                 // use a prefix string (separated by colons)
                                 // for each line of output; if it calls
                                 // another function, that may also use its
                                 // prefix which is then printed after the one
                                 // of the calling function. Since output from
                                 // functions which are nested deep below is
                                 // usually not as important as top-level
                                 // output, you can give the deallog variable
                                 // a maximal depth of nested output for
                                 // output to console and file. The depth zero
                                 // which we gave here means that no output is
                                 // written. By changing it you can get more
                                 // information about the innards of the
                                 // library.
int main (int argc, char *argv[])
{
  deallog.depth_console (0);

  clock_t     start, end;

  start = clock();

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);
  {
    Step4<3> laplace_problem_3d;
    laplace_problem_3d.run ();
  }

  end = clock();
  cout<< "%%%%%% Rechenzeit overall = " << (double)(end-start)/CLOCKS_PER_SEC <<std::endl;

  return 0;
}
