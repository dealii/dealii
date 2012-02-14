/* Authors: Joerg Frohne, Texas A&M University and                */
/*                        University of Siegen, 2011, 2012        */
/*          Wolfgang Bangerth, Texas A&M University, 2012         */

/*    $Id: step-4.cc 24093 2011-08-16 13:58:12Z bangerth $        */
/*                                                                */
/*    Copyright (C) 2011, 2012 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

                                 // @sect3{Include files}

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vectors.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <list>


namespace Step41
{
  using namespace dealii;

				   // @sect3{The <code>Step41</code> class template}

				   // This class supplies all function and
				   // variables to an obstacle problem. The
				   // update_solution_and_constraints function and the
				   // ConstaintMatrix are important for the
				   // handling of the active set as we see
				   // later.

  template <int dim>
  class ObstacleProblem
  {
    public:
      ObstacleProblem ();
      void run ();

    private:
      void make_grid ();
      void setup_system();
      void assemble_system ();
      void assemble_mass_matrix (TrilinosWrappers::SparseMatrix &mass_matrix);
      void update_solution_and_constraints ();
      void solve ();
      void output_results (const unsigned int iteration) const;

      Triangulation<dim>   triangulation;
      FE_Q<dim>            fe;
      DoFHandler<dim>      dof_handler;
      ConstraintMatrix     constraints;
      IndexSet             active_set;

      TrilinosWrappers::SparseMatrix system_matrix;
      TrilinosWrappers::SparseMatrix complete_system_matrix;

      TrilinosWrappers::Vector       solution;
      TrilinosWrappers::Vector       system_rhs;
      TrilinosWrappers::Vector       complete_system_rhs;
      TrilinosWrappers::Vector       diagonal_of_mass_matrix;
  };


				   // @sect3{Right hand side and boundary values}

  template <int dim>
  class RightHandSide : public Function<dim>
  {
    public:
      RightHandSide () : Function<dim>() {}

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;
  };

  template <int dim>
  class BoundaryValues : public Function<dim>
  {
    public:
      BoundaryValues () : Function<dim>() {}

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;
  };

  template <int dim>
  class Obstacle : public Function<dim>
  {
    public:
      Obstacle () : Function<dim>() {}

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;
  };



				   // For this example, we choose as right hand
				   // side function a constant force density
				   // like the gravitation attraction.
  template <int dim>
  double RightHandSide<dim>::value (const Point<dim> &p,
				    const unsigned int component) const
  {
    Assert (component == 0, ExcNotImplemented());

    return -10;
  }


				   // As boundary values, we choose the zero.
  template <int dim>
  double BoundaryValues<dim>::value (const Point<dim> &p,
				     const unsigned int component) const
  {
    Assert (component == 0, ExcNotImplemented());

    return 0;
  }


				   // The obstacle function describes a cascaded
				   // barrier. So if the gravitation attraction
				   // pulls the membrane down it blows over the
				   // steps.
  template <int dim>
  double Obstacle<dim>::value (const Point<dim> &p,
			       const unsigned int component) const
  {
    Assert (component == 0, ExcNotImplemented());

    if (p (0) < -0.5)
      return -0.2;
    else if (p (0) >= -0.5 && p (0) < 0.0)
      return -0.4;
    else if (p (0) >= 0.0 && p (0) < 0.5)
      return -0.6;
    else
      return -0.8;
  }



				   // @sect3{Implementation of the <code>ObstacleProblem</code> class}


				   // @sect4{ObstacleProblem::ObstacleProblem}

  template <int dim>
  ObstacleProblem<dim>::ObstacleProblem ()
		  :
		  fe (1),
		  dof_handler (triangulation)
  {}


				   // @sect4{ObstacleProblem::make_grid}

				   // We solve our obstacle problem on the square
				   // $[-1,1]\times [-1,1]$ in 2D.
  template <int dim>
  void ObstacleProblem<dim>::make_grid ()
  {
    GridGenerator::hyper_cube (triangulation, -1, 1);
    triangulation.refine_global (7);

    std::cout << "Number of active cells: "
	      << triangulation.n_active_cells()
	      << std::endl
	      << "Total number of cells: "
	      << triangulation.n_cells()
    	      << std::endl;
  }

				   // @sect4{ObstacleProblem::setup_system}

  template <int dim>
  void ObstacleProblem<dim>::setup_system ()
  {
    dof_handler.distribute_dofs (fe);
    active_set.set_size (dof_handler.n_dofs());

    std::cout << "Number of degrees of freedom: "
	      << dof_handler.n_dofs()
	      << std::endl
	      << std::endl;

    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      BoundaryValues<dim>(),
					      constraints);
    constraints.close ();

    CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler,
				     c_sparsity,
				     constraints,
				     false);

    system_matrix.reinit (c_sparsity);
    complete_system_matrix.reinit (c_sparsity);

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
    complete_system_rhs.reinit (dof_handler.n_dofs());

				     // to compute the factor which is used
				     // to scale the residual. You can consider
				     // this diagonal matrix as the discretization
				     // of a lagrange multiplier for the
				     // contact force
    TrilinosWrappers::SparseMatrix mass_matrix;
    mass_matrix.reinit (c_sparsity);
    assemble_mass_matrix (mass_matrix);
    diagonal_of_mass_matrix.reinit (dof_handler.n_dofs());
    for (unsigned int j=0; j<solution.size (); j++)
      diagonal_of_mass_matrix (j) = mass_matrix.diag_element (j);
  }


				   // @sect4{ObstacleProblem::assemble_system}


				   // At once with assembling the system matrix and
				   // right-hand-side we apply the constraints
				   // to our system. The constraint consists not
				   // only of the zero Dirichlet boundary values,
				   // in addition they contain the obstacle values.
				   // The update_solution_and_constraints function are used
				   // to fill the ConstraintMatrix.
  template <int dim>
  void ObstacleProblem<dim>::assemble_system ()
  {
    std::cout << "   Assembling system..." << std::endl;

    system_matrix = 0;
    system_rhs    = 0;

    const QGauss<dim>         quadrature_formula(2);
    const RightHandSide<dim>  right_hand_side;

    FEValues<dim>             fe_values (fe, quadrature_formula,
					update_values   | update_gradients |
					update_quadrature_points |
					update_JxW_values);

    const unsigned int        dofs_per_cell = fe.dofs_per_cell;
    const unsigned int        n_q_points    = quadrature_formula.size();

    FullMatrix<double>        cell_matrix (dofs_per_cell, dofs_per_cell);
    TrilinosWrappers::Vector  cell_rhs (dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      {
	fe_values.reinit (cell);
	cell_matrix = 0;
	cell_rhs = 0;

	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
				     fe_values.shape_grad (j, q_point) *
				     fe_values.JxW (q_point));

	      cell_rhs(i) += (fe_values.shape_value (i, q_point) *
			      right_hand_side.value (fe_values.quadrature_point (q_point)) *
			      fe_values.JxW (q_point));
	    }

	cell->get_dof_indices (local_dof_indices);

					 // This function apply the constraints
					 // to the system matrix and system rhs.
					 // The true parameter is set to make sure
					 // that the system rhs contains correct
					 // values in the rows with inhomogeneity
					 // constraints.
	constraints.distribute_local_to_global (cell_matrix,
						cell_rhs,
						local_dof_indices,
						system_matrix,
						system_rhs,
						true);
      }
  }


  template <int dim>
  void
  ObstacleProblem<dim>::
  assemble_mass_matrix (TrilinosWrappers::SparseMatrix &mass_matrix)
  {
    const QTrapez<dim>        quadrature_formula;
    FEValues<dim>             fe_values (fe,
					 quadrature_formula,
					 update_values   |
					 update_quadrature_points |
					 update_JxW_values);

    const unsigned int        dofs_per_cell = fe.dofs_per_cell;
    const unsigned int        n_q_points    = quadrature_formula.size();

    FullMatrix<double>        cell_matrix (dofs_per_cell, dofs_per_cell);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      {
	fe_values.reinit (cell);
	cell_matrix = 0;

	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (fe_values.shape_value (i, q_point) *
				   fe_values.shape_value (j, q_point) *
				   fe_values.JxW (q_point));

	cell->get_dof_indices (local_dof_indices);

					 // This function apply the constraints
					 // to the system matrix and system rhs.
					 // The true parameter is set to make sure
					 // that the system rhs contains correct
					 // values in the rows with inhomogeneity
					 // constraints.
	constraints.distribute_local_to_global (cell_matrix,
						local_dof_indices,
						mass_matrix);
      }
  }

				   // @sect4{ObstacleProblem::update_solution_and_constraints}

				   // Updating of the active set which means to
				   // set a inhomogeneity constraint in the
				   // ConstraintMatrix. At the same time we set
				   // the solution to the correct value - the obstacle value.
				   // To control the active set we use the vector
				   // active_set which contains a zero in a component
				   // that is not in the active set and elsewise a
				   // one. With the output file you can visualize it.
  template <int dim>
  void
  ObstacleProblem<dim>::update_solution_and_constraints ()
  {
    std::cout << "   Updating active set..." << std::endl;

    const Obstacle<dim>     obstacle;
    unsigned int            counter_contact_constraints = 0;


    TrilinosWrappers::Vector       force_residual (dof_handler.n_dofs());
    complete_system_matrix.residual (force_residual,
				     solution, complete_system_rhs);
    force_residual *= -1;

    constraints.clear();

				     // to find and supply the constraints for the
				     // obstacle condition
    active_set.clear ();
    const double c = 100.0;

    std::vector<bool> dof_touched (dof_handler.n_dofs(),
				   false);

    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
	{
	  Assert (dof_handler.get_fe().dofs_per_cell ==
		  GeometryInfo<dim>::vertices_per_cell,
		  ExcNotImplemented());

	  const unsigned int dof_index = cell->vertex_dof_index (v,0);

	  if (dof_touched[dof_index] == true)
	    continue;

					   // the local row where
	  const double obstacle_value = obstacle.value (cell->vertex(v));
	  const double solution_value = solution (dof_index);

					   // To decide which dof belongs to the
					   // active-set. For that we scale the
					   // residual-vector with the cell-size and
					   // the diag-entry of the mass-matrix.

					   // TODO: I have to check the condition

	  if (force_residual (dof_index) +
	      c * diagonal_of_mass_matrix(dof_index) * (obstacle_value - solution_value)
	      >
	      0)
	    {
	      active_set.add_index (dof_index);
	      constraints.add_line (dof_index);
	      constraints.set_inhomogeneity (dof_index, obstacle_value);

	      solution (dof_index) = obstacle_value;
					       // the residual of the non-contact
					       // part of the system serves as an
					       // additional control which is not
					       // necessary for for the primal-dual
					       // active set strategy
	      force_residual (dof_index) = 0;

	      dof_touched[dof_index] = true;
	    }
	}
    std::cout << "      Size of active set: " << active_set.n_elements()
	      << std::endl;

				     // To supply the boundary values of the
				     // dirichlet-boundary in constraints
    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      BoundaryValues<dim>(),
					      constraints);
    constraints.close ();

    std::cout << "   Residual of the non-contact part of the system: "
	      << force_residual.l2_norm()
	      << std::endl;

  }

				   // @sect4{ObstacleProblem::solve}

  template <int dim>
  void ObstacleProblem<dim>::solve ()
  {
    std::cout << "   Solving system..." << std::endl;

    ReductionControl        reduction_control (100, 1e-12, 1e-3);
    SolverCG<TrilinosWrappers::Vector>   solver (reduction_control);
    TrilinosWrappers::PreconditionAMG precondition;
    precondition.initialize (system_matrix);

    solver.solve (system_matrix, solution, system_rhs, precondition);
    constraints.distribute (solution);

    std::cout << "      Error: " << reduction_control.initial_value()
	      << " -> " << reduction_control.last_value()
	      << " in "
	      <<  reduction_control.last_step()
	      << " CG iterations."
	      << std::endl;
  }


				   // @sect4{ObstacleProblem::output_results}

				   // We use the vtk-format for the
				   // output.  The file contains the
				   // displacement and a numerical
				   // represenation of the active
				   // set. The function looks standard
				   // but note that we can add an
				   // IndexSet object to the DataOut
				   // object in exactly the same way
				   // as a regular solution vector: it
				   // is simply interpreted as a
				   // function that is either zero
				   // (when a degree of freedom is not
				   // part of the IndexSet) or one (if
				   // it is).
  template <int dim>
  void ObstacleProblem<dim>::output_results (const unsigned int iteration) const
  {
    std::cout << "   Writing graphical output..." << std::endl;

    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "displacement");
    data_out.add_data_vector (active_set, "active_set");

    data_out.build_patches ();

    std::ofstream output_vtk ((std::string("output_") +
			       Utilities::int_to_string (iteration, 3) +
			       ".vtk").c_str ());
    data_out.write_vtk (output_vtk);
  }



				   // @sect4{ObstacleProblem::run}

				   // This is the function which has
				   // the top-level control over
				   // everything.  It is not very
				   // long, and in fact rather
				   // straightforward: in every
				   // iteration of the active set
				   // method, we assemble the linear
				   // system, solve it, update the
				   // active set and project the
				   // solution back to the feasible
				   // set, and then output the
				   // results. The iteration is
				   // terminated whenever the active
				   // set has not changed in the
				   // previous iteration.
				   //
				   // The only trickier part is that
				   // we have to save the linear
				   // system (i.e., the matrix and
				   // right hand side) after
				   // assembling it in the first
				   // iteration. The reason is that
				   // this is the only step where we
				   // can access the linear system as
				   // built without any of the contact
				   // constraints active. We need this
				   // to compute the residual of the
				   // solution at other iterations,
				   // but in other iterations that
				   // linear system we form has the
				   // rows and columns that correspond
				   // to constrained degrees of
				   // freedom eliminated, and so we
				   // can no longer access the full
				   // residual of the original
				   // equation.
  template <int dim>
  void ObstacleProblem<dim>::run ()
  {
    make_grid();
    setup_system ();

    IndexSet active_set_old (active_set);
    for (unsigned int iteration=0; iteration<=solution.size (); ++iteration)
      {
	std::cout << "Newton iteration " << iteration << std::endl;

	assemble_system ();

	if (iteration == 0)
	  {
	    complete_system_matrix.copy_from (system_matrix);
	    complete_system_rhs = system_rhs;
	  }

	solve ();
	update_solution_and_constraints ();
	output_results (iteration);

	if (active_set == active_set_old)
	  break;

	active_set_old = active_set;

	std::cout << std::endl;
      }
  }
}


                                 // @sect3{The <code>main</code> function}

				 // And this is the main function. It
				 // follows the pattern of all other
				 // main functions. The call to
				 // initialize MPI exists because the
				 // Trilinos library upon which we
				 // build our linear solvers in this
				 // program requires it.
int main (int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step41;

      deallog.depth_console (0);

      Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

      ObstacleProblem<2> obstacle_problem;
      obstacle_problem.run ();
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
