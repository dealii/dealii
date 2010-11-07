/* $Id$ */
/* Author: Wolfgang Bangerth, Texas A&M University, 2009, 2010 */
/*         Timo Heister, University of Goettingen, 2009, 2010 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2009, 2010 by Timo Heister and the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


                                 // @sect3{Include files}
				 //
				 // Most of the include files we need for this
				 // program have already been discussed in
				 // previous programs. In particular, all of
				 // the following should already be familiar
				 // friends:
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/solver_cg.h>
#include <lac/constraint_matrix.h>
#include <lac/compressed_simple_sparsity_pattern.h>

#include <lac/petsc_parallel_sparse_matrix.h>
#include <lac/petsc_parallel_vector.h>
#include <lac/petsc_solver.h>
#include <lac/petsc_precondition.h>

#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <fe/fe_q.h>
#include <numerics/vectors.h>
#include <numerics/data_out.h>
#include <numerics/error_estimator.h>

				 // The following, however, will be new or be
				 // used in new roles. Let's walk through
				 // them. The first of these will provide the
				 // tools of the Utilities::System namespace
				 // that we will use to query things like the
				 // number of processors associated with the
				 // current MPI universe, or the number within
				 // this universe the processor this job runs
				 // on is:
#include <base/utilities.h>
				 // The next one provides a class,
				 // ConditionOStream that allows us to write
				 // code that would output things to a stream
				 // (such as <code>std::cout</code> on every
				 // processor but throws the text away on all
				 // but one of them. We could achieve the same
				 // by simply putting an <code>if</code>
				 // statement in front of each place where we
				 // may generate output, but this doesn't make
				 // the code any prettier. In addition, the
				 // condition whether this processor should or
				 // should not produce output to the screen is
				 // the same every time -- and consequently it
				 // should be simple enough to put it into the
				 // statements that generate output itself.
#include <base/conditional_ostream.h>
				 // After these preliminaries, here is where
				 // it becomes more interesting. As mentioned
				 // in the @ref distributed module, one of the
				 // fundamental truths of solving problems on
				 // large numbers of processors is that there
				 // is no way for any processor to store
				 // everything (e.g. information about all
				 // cells in the mesh, all degrees of freedom,
				 // or the values of all elements of the
				 // solution vector). Rather, every processor
				 // will <i>own</i> a few of each of these
				 // and, if necessary, may <i>know</i> about a
				 // few more, for example the ones that are
				 // located on cells adjacent to the ones this
				 // processor owns itself. We typically call
				 // the latter <i>ghost cells</i>, <i>ghost
				 // nodes</i> or <i>ghost elements of a
				 // vector</i>. The point of this discussion
				 // here is that we need to have a way to
				 // indicate which elements a particular
				 // processor owns or need to know of. This is
				 // the realm of the IndexSet class: if there
				 // are a total of $N$ cells, degrees of
				 // freedom, or vector elements, associated
				 // with (non-negative) integral indices
				 // $[0,N)$, then both the set of elements the
				 // current processor owns as well as the
				 // (possibly larger) set of indices it needs
				 // to know about are subsets of the set
				 // $[0,N)$. IndexSet is a class that stores
				 // subsets of this set in an efficient
				 // format:
#include <base/index_set.h>
				 // The next header file is necessary for a
				 // single function,
				 // SparsityTools::distribute_sparsity_pattern. The
				 // role of this function will be explained
				 // below.
#include <lac/sparsity_tools.h>
				 // The final two, new header files provide
				 // the class
				 // parallel::distributed::Triangulation that
				 // provides meshes distributed across a
				 // potentially very large number of
				 // processors, while the second provides the
				 // namespace
				 // parallel::distributed::GridRefinement that
				 // offers functions that can adaptively
				 // refine such distributed meshes:
#include <distributed/tria.h>
#include <distributed/grid_refinement.h>

#include <fstream>
#include <iostream>

using namespace dealii;

                                 // @sect3{The <code>LaplaceProblem</code> class template}

				 // Next let's declare the main class of this
				 // program. Its structure is almost exactly
				 // that of the step-6 tutorial program. The
				 // only significant differences are:
				 // - The <code>mpi_communicator</code>
				 //   variable that describes the set of
				 //   processors we want this code to run
				 //   on. In practice, this will be
				 //   MPI_COMM_WORLD, i.e. all processors the
				 //   batch scheduling system has assigned to
				 //   this particular job.
				 // - The presence of the <code>pcout</code>
				 //   variable of type ConditionOStream.
				 // - The obvious use of
				 //   parallel::distributed::Triangulation
				 //   instead of Triangulation.
				 // - The fact that all matrices and vectors
				 //   are now distributed. We use their
				 //   PETScWrapper versions for this since
				 //   deal.II's own classes do not provide
				 //   %parallel functionality.
				 // - The presence of two IndexSet objects
				 //   that denote which sets of degrees of
				 //   freedom (and associated elements of
				 //   solution and right hand side vectors) we
				 //   own on the current processor and which
				 //   we need (as ghost elements) for the
				 //   algorithms in this program to work.
template <int dim>
class LaplaceProblem
{
  public:
    LaplaceProblem ();
    ~LaplaceProblem ();

    void run (const unsigned int initial_global_refine);

  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    MPI_Comm mpi_communicator;

    parallel::distributed::Triangulation<dim>   triangulation;

    DoFHandler<dim>      dof_handler;
    FE_Q<dim>            fe;

    IndexSet             locally_owned_dofs;
    IndexSet             locally_relevant_dofs;

    ConstraintMatrix     constraints;

    PETScWrappers::MPI::SparseMatrix system_matrix;
    PETScWrappers::MPI::Vector locally_relevant_solution;
    PETScWrappers::MPI::Vector system_rhs;

    ConditionalOStream                pcout;
};


                                 // @sect3{The <code>LaplaceProblem</code> class implementation}

                                 // @sect4{Constructors and destructors}

				 // Constructors and destructors are rather
				 // trivial. In addition to what we do in
				 // step-6, we set the set of processors we
				 // want to work on to all machines available
				 // (MPI_COMM_WORLD); ask the triangulation to
				 // ensure that the mesh remains smooth and
				 // free to refined islands, for example; and
				 // initialize the <code>pcout</code> variable
				 // to only allow processor zero to output
				 // anything:
template <int dim>
LaplaceProblem<dim>::LaplaceProblem ()
		:
		mpi_communicator (MPI_COMM_WORLD),
		triangulation (mpi_communicator,
			       typename Triangulation<dim>::MeshSmoothing
			       (Triangulation<dim>::smoothing_on_refinement |
				Triangulation<dim>::smoothing_on_coarsening)),
		dof_handler (triangulation),
                fe (2),
		pcout (std::cout,
		       (Utilities::System::
			get_this_mpi_process(mpi_communicator)
			== 0))
{}



template <int dim>
LaplaceProblem<dim>::~LaplaceProblem ()
{
  dof_handler.clear ();
}


                                 // @sect4{LaplaceProblem::setup_system}

				 // The following function is, arguably, the
				 // most interesting one in the entire program
				 // since it goes to the heart of what
				 // distinguishes %parallel step-40 from
				 // sequential step-6.
				 //
				 // At the top we do what we always do: tell
				 // the DoFHandler object to distribute
				 // degrees of freedom. Since the
				 // triangulation we use here is distributed,
				 // the DoFHandler object is smart enough to
				 // recognize that on each processor it can
				 // only distribute degrees of freedom on
				 // cells it owns; this is followed by an
				 // exchange step in which processors tell
				 // each other about degrees of freedom on
				 // ghost cell. The result is a DoFHandler
				 // that knows about the degrees of freedom on
				 // locally owned cells and ghost cells
				 // (i.e. cells adjacent to locally owned
				 // cells) but nothing about cells that are
				 // further away, consistent with the basic
				 // philosophy of distributed computing that
				 // no processor can know everything.
template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  locally_owned_dofs = dof_handler.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (dof_handler,
					   locally_relevant_dofs);

				   // note in class doc!
  locally_relevant_solution.reinit (mpi_communicator,
				    locally_owned_dofs,
				    locally_relevant_dofs);
  system_rhs.reinit (mpi_communicator,
		     dof_handler.n_dofs(),
		     dof_handler.n_locally_owned_dofs());

  locally_relevant_solution = 0;
  system_rhs = 0;

  constraints.clear ();
  constraints.reinit (locally_relevant_dofs);

  DoFTools::make_hanging_node_constraints (dof_handler, constraints);
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    constraints);
  constraints.close ();


  CompressedSimpleSparsityPattern csp (dof_handler.n_dofs(),
				       dof_handler.n_dofs(),
				       locally_relevant_dofs);
  DoFTools::make_sparsity_pattern (dof_handler,
				   csp,
				   constraints, false);
  SparsityTools::distribute_sparsity_pattern (csp,
					      dof_handler.n_locally_owned_dofs_per_processor(),
					      mpi_communicator,
					      locally_relevant_dofs);
  system_matrix.reinit (mpi_communicator,
			csp,
			dof_handler.n_locally_owned_dofs_per_processor(),
			dof_handler.n_locally_owned_dofs_per_processor(),
			Utilities::System::get_this_mpi_process(mpi_communicator));
}



template <int dim>
void LaplaceProblem<dim>::assemble_system ()
{
  const QGauss<dim>  quadrature_formula(3);

  FEValues<dim> fe_values (fe, quadrature_formula,
			   update_values    |  update_gradients |
			   update_quadrature_points |
			   update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    if (cell->subdomain_id() == triangulation.locally_owned_subdomain())
      {
	cell_matrix = 0;
	cell_rhs = 0;

	fe_values.reinit (cell);

	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
				     fe_values.shape_grad(j,q_point) *
				     fe_values.JxW(q_point));

	      cell_rhs(i) += (fe_values.shape_value(i,q_point) *
			      (fe_values.quadrature_point(q_point)[1]
			       >
			       0.5+0.25*sin(4.0*numbers::PI*fe_values.quadrature_point(q_point)[0])
			       ? 1 : -1) *
			      fe_values.JxW(q_point));
	    }

	cell->get_dof_indices (local_dof_indices);
	constraints.distribute_local_to_global (cell_matrix,
						cell_rhs,
						local_dof_indices,
						system_matrix,
						system_rhs);
      }

  system_matrix.compress ();
  system_rhs.compress ();
}





template <int dim>
void LaplaceProblem<dim>::solve ()
{
  PETScWrappers::MPI::Vector
    completely_distributed_solution (mpi_communicator,
				     dof_handler.n_dofs(),
				     dof_handler.n_locally_owned_dofs());

  SolverControl solver_control (dof_handler.n_dofs(), 1e-12);

  PETScWrappers::SolverCG solver(solver_control, mpi_communicator);
  PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);

  solver.solve (system_matrix, completely_distributed_solution, system_rhs,
		preconditioner);

  pcout << "  Solved in " << solver_control.last_step()
	<< " iterations." << std::endl;

  constraints.distribute (completely_distributed_solution);

  locally_relevant_solution = completely_distributed_solution;
}



template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
  KellyErrorEstimator<dim>::estimate (dof_handler,
				      QGauss<dim-1>(3),
				      typename FunctionMap<dim>::type(),
				      locally_relevant_solution,
 				      estimated_error_per_cell);
  parallel::distributed::GridRefinement::
    refine_and_coarsen_fixed_number (triangulation,
				     estimated_error_per_cell,
				     0.3, 0.03);
  triangulation.execute_coarsening_and_refinement ();
}




template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (locally_relevant_solution, "u");

  Vector<float> subdomain (triangulation.n_active_cells());
				   // could just fill entire vector with subdomain_id()
  {
    unsigned int index = 0;
    for (typename Triangulation<dim>::active_cell_iterator
	   cell = triangulation.begin_active();
	 cell != triangulation.end(); ++cell, ++index)
      subdomain(index) = (cell->is_ghost() || cell->is_artificial()
			  ?
			  -1
			  :
			  cell->subdomain_id());
  }
  data_out.add_data_vector (subdomain, "subdomain");
  data_out.build_patches ();

  const std::string filename = ("solution-" +
				Utilities::int_to_string (cycle, 2) +
				"." +
				Utilities::int_to_string
				(triangulation.locally_owned_subdomain(), 4));

  std::ofstream output ((filename + ".vtu").c_str());
  data_out.write_vtu (output);

  if (Utilities::System::get_this_mpi_process(mpi_communicator) == 0)
    {
      std::vector<std::string> filenames;
      for (unsigned int i=0;
	   i<Utilities::System::get_n_mpi_processes(mpi_communicator);
	   ++i)
	filenames.push_back ("solution-" +
			     Utilities::int_to_string (cycle, 2) +
			     "." +
			     Utilities::int_to_string (i, 4) +
			     ".vtu");

      std::ofstream master_output ((filename + ".pvtu").c_str());
      data_out.write_pvtu_record (master_output, filenames);
    }
}



template <int dim>
void LaplaceProblem<dim>::run (const unsigned int initial_global_refine)
{
  const unsigned int n_cycles = 12;
  for (unsigned int cycle=0; cycle<n_cycles; ++cycle)
    {
      pcout << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
	{
	  GridGenerator::hyper_cube (triangulation);
	  triangulation.refine_global (initial_global_refine);
	}
      else
	refine_grid ();

      setup_system ();

      pcout << "   Number of active cells:       "
	    << triangulation.n_global_active_cells()
	    << std::endl
	    << "   Number of degrees of freedom: "
	    << dof_handler.n_dofs()
	    << std::endl;
      
      assemble_system ();
      solve ();

      if (Utilities::System::get_n_mpi_processes(mpi_communicator) <= 100)
	output_results (cycle);

      pcout << std::endl;
    }
}



int main(int argc, char *argv[])
{
  try
    {
      PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
      deallog.depth_console (0);

      int refine=5;
      if (argc>1)
	{
	  refine = (unsigned int)Utilities::string_to_int(argv[1]);
	}

      {
	LaplaceProblem<2> laplace_problem_2d;
	laplace_problem_2d.run (refine);
      }

      PetscFinalize();
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
