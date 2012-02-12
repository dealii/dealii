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
				   // projection_active_set function and the
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
      void projection_active_set ();
      void solve ();
      void output_results (const unsigned int iteration) const;

      Triangulation<dim>   triangulation;
      FE_Q<dim>            fe;
      DoFHandler<dim>      dof_handler;
      ConstraintMatrix     constraints;
      IndexSet             active_set;

      TrilinosWrappers::SparseMatrix system_matrix;
      TrilinosWrappers::SparseMatrix system_matrix_complete;

      TrilinosWrappers::Vector       solution;
      TrilinosWrappers::Vector       system_rhs;
      TrilinosWrappers::Vector       system_rhs_complete;
      TrilinosWrappers::Vector       force_residual;
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
    system_matrix_complete.reinit (c_sparsity);

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
    system_rhs_complete.reinit (dof_handler.n_dofs());
    force_residual.reinit (dof_handler.n_dofs());

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
				   // The projection_active_set function are used
				   // to fill the ConstraintMatrix.
  template <int dim>
  void ObstacleProblem<dim>::assemble_system ()
  {
    std::cout << "   Assembling system..." << std::endl;

    QGauss<dim>  quadrature_formula(2);

    const RightHandSide<dim> right_hand_side;

    FEValues<dim> fe_values (fe, quadrature_formula,
			     update_values   | update_gradients |
			     update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    TrilinosWrappers::Vector       cell_rhs (dofs_per_cell);

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
	constraints.distribute_local_to_global (cell_matrix, cell_rhs,
						local_dof_indices,
						system_matrix, system_rhs, true);
      }
  }

  template <int dim>
  void ObstacleProblem<dim>::assemble_mass_matrix (TrilinosWrappers::SparseMatrix &mass_matrix)
  {
    QTrapez<dim>  quadrature_formula;

    FEValues<dim> fe_values (fe, quadrature_formula,
			     update_values   | update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

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

				   // @sect4{ObstacleProblem::projection_active_set}

				   // Updating of the active set which means to
				   // set a inhomogeneity constraint in the
				   // ConstraintMatrix. At the same time we set
				   // the solution to the correct value - the obstacle value.
				   // To control the active set we use the vector
				   // active_set which contains a zero in a component
				   // that is not in the active set and elsewise a
				   // one. With the output file you can visualize it.
  template <int dim>
  void ObstacleProblem<dim>::projection_active_set ()
  {
    std::cout << "   Updating active set..." << std::endl;

    const Obstacle<dim>     obstacle;
    std::vector<bool>       vertex_touched (triangulation.n_vertices(),
					    false);
    unsigned int            counter_contact_constraints = 0;

    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

    constraints.clear();

				     // to find and supply the constraints for the
				     // obstacle condition
    active_set.clear ();
    const double c = 100.0;
    for (; cell!=endc; ++cell)
      for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
	{
	  unsigned int index_x = cell->vertex_dof_index (v,0);

					   // the local row where
	  Point<dim> point (cell->vertex (v)[0], cell->vertex (v)[1]);
	  double obstacle_value = obstacle.value (point);
	  double solution_index_x = solution (index_x);

					   // To decide which dof belongs to the
					   // active-set. For that we scale the
					   // residual-vector with the cell-size and
					   // the diag-entry of the mass-matrix.

					   // TODO: I have to check the condition
	  if (force_residual (index_x) +
	      diagonal_of_mass_matrix (index_x)*c*(obstacle_value - solution_index_x) > 0)
	    {
	      constraints.add_line (index_x);
	      constraints.set_inhomogeneity (index_x, obstacle_value);
	      solution (index_x) = obstacle_value;
	      active_set.add_index (index_x);

	      if (vertex_touched[cell->vertex_index(v)] == false)
		{
		  vertex_touched[cell->vertex_index(v)] = true;
		  counter_contact_constraints += 1;
		}
	    }
	}
    std::cout << "      Size of active set: " << counter_contact_constraints
	      << std::endl;

				     // To supply the boundary values of the
				     // dirichlet-boundary in constraints
    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      BoundaryValues<dim>(),
					      constraints);
    constraints.close ();
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

				   // We use the vtk-format for the output.
				   // The file contains the displacement,
				   // the residual and active set vectors.
  template <int dim>
  void ObstacleProblem<dim>::output_results (const unsigned int iteration) const
  {
    std::cout << "   Writing graphical output..." << std::endl;

    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "displacement");
    data_out.add_data_vector (force_residual, "residual");

    Vector<double> numerical_active_set (dof_handler.n_dofs());
    active_set.fill_binary_vector (numerical_active_set);
    data_out.add_data_vector (numerical_active_set, "active_set");

    data_out.build_patches ();

    std::ofstream output_vtk ((std::string("output_") +
			       Utilities::int_to_string (iteration) +
			       ".vtk").c_str ());
    data_out.write_vtk (output_vtk);
  }



				   // @sect4{ObstacleProblem::run}

				   // This is the function which has the
				   // top-level control over everything.
				   // Here the active set method is implemented.

				   // TODO: I have to compare it with the algorithm
				   // in the Wohlmuth-paper
  template <int dim>
  void ObstacleProblem<dim>::run ()
  {
    make_grid();
    setup_system ();

				     // TODO: can't some of this be
				     // merged with the first Newton
				     // iteration?
    std::cout << "Initial start-up step" << std::endl;

    ConstraintMatrix constraints_complete (constraints);
    assemble_system ();
    solve ();

				     // to save the system_matrix and the
				     // rhs to compute the residual in every
				     // step of the active-set-iteration
    system_matrix_complete.copy_from (system_matrix);
    system_rhs_complete = system_rhs;

				     //TODO: use system_matrix_complete.residual
    force_residual = 0;
    force_residual -= system_rhs_complete;
    system_matrix_complete.vmult_add  (force_residual, solution);

				     // to compute a start active set
    projection_active_set ();

    std::cout << std::endl;

    IndexSet active_set_old (active_set);
    for (unsigned int iteration=1; iteration<=solution.size (); ++iteration)
      {
	std::cout << "Newton iteration " << iteration << std::endl;

	system_matrix = 0;
	system_rhs = 0;

	assemble_system ();
	solve ();

				     //TODO: use system_matrix_complete.residual
	force_residual = 0;
	force_residual -= system_rhs_complete;
	system_matrix_complete.vmult_add  (force_residual, solution);

	projection_active_set ();

	for (unsigned int k = 0; k<solution.size (); k++)
	  if (active_set.is_element (k))
	    force_residual (k) = 0;

	output_results (iteration);

					 // the residual of the non-contact part
					 // of the system serves as an additional
					 // control which is not necassary for
					 // for the primal-dual active set strategy
	std::cout << "   Residual of the non-contact part of the system: "
		  << force_residual.l2_norm()
		  << std::endl;

					 // if both the old and the new
					 // active set are identical the
					 // computation stops
	if (active_set == active_set_old)
	  break;
	active_set_old = active_set;

	std::cout << std::endl;
      }
  }
}


                                 // @sect3{The <code>main</code> function}

				 // And this is the main function. It also
				 // looks mostly like in step-3, but if you
				 // look at the code below, note how we first
				 // create a variable of type
				 // <code>ObstacleProblem@<2@></code> (forcing
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
