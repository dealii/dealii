/* $Id$ */
/* Author: Wolfgang Bangerth, Texas A&M University, 2009 */
/*         Timo Heister, University of Goettingen, 2009 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2009, 2010 by Timo Heister and the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

// comment to use PETSc:
//#define USE_TRILINOS


#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <base/utilities.h>
#include <base/conditional_ostream.h>
#include <base/index_set.h>
#include <base/timer.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/constraint_matrix.h>
#include <lac/compressed_simple_sparsity_pattern.h>
#include <lac/sparsity_tools.h>

#ifdef USE_TRILINOS
#include <lac/trilinos_precondition.h>
#include <lac/trilinos_sparse_matrix.h>
#include <lac/trilinos_sparsity_pattern.h>
#include <lac/trilinos_vector.h>
#else
#include <lac/petsc_parallel_sparse_matrix.h>
#include <lac/petsc_parallel_vector.h>
#include <lac/petsc_solver.h>
#include <lac/petsc_precondition.h>
#endif

#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_out.h>
#include <grid/filtered_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <fe/fe_q.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/error_estimator.h>

#include <distributed/tria.h>
#include <distributed/grid_refinement.h>
#include <distributed/solution_transfer.h>

#include <fstream>
#include <iostream>

#include "timeblock.h"

using namespace dealii;

void print_it(Utilities::System::MinMaxAvg & result)
{
  std::cout// << "sum: " << result.sum
    << " avg: " << (long)result.avg/1024
	  << " min: " << (long)result.min/1024 << " @" << result.min_index
	  << " max: " << (long)result.max/1024 << " @" << result.max_index
	  << std::endl;
}

void print_memory_stats()
{
  int myid = Utilities::System::get_this_mpi_process(MPI_COMM_WORLD);
  Utilities::System::MemoryStats stats;
  Utilities::System::get_memory_stats(stats);
  Utilities::System::MinMaxAvg r;
  Utilities::System::calculate_collective_mpi_min_max_avg(MPI_COMM_WORLD, stats.VmPeak, r);
  if (myid==0)
    {
      std::cout << "MEM: VmPeak: ";
      print_it(r);
    }
  Utilities::System::calculate_collective_mpi_min_max_avg(MPI_COMM_WORLD, stats.VmSize, r);
  if (myid==0)
    {
      std::cout << "MEM: VmSize: ";
      print_it(r);
    }
  Utilities::System::calculate_collective_mpi_min_max_avg(MPI_COMM_WORLD, stats.VmHWM, r);
  if (myid==0)
    {
      std::cout << "MEM: VmHWM:  ";
      print_it(r);
    }
  Utilities::System::calculate_collective_mpi_min_max_avg(MPI_COMM_WORLD, stats.VmRSS, r);
  if (myid==0)
    {
      std::cout << "MEM: VmRSS:  ";
      print_it(r);
    }
}


template <int dim>
class LaplaceProblem
{
  public:
    LaplaceProblem ();
    ~LaplaceProblem ();

    void run (unsigned int refine);

  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    parallel::distributed::Triangulation<dim>   triangulation;

    DoFHandler<dim>      dof_handler;
    FE_Q<dim>            fe;

    IndexSet             locally_owned_dofs;
    IndexSet             locally_relevant_dofs;

    ConstraintMatrix     constraints;

#ifdef USE_TRILINOS
    TrilinosWrappers::SparseMatrix    system_matrix;
    TrilinosWrappers::MPI::Vector     solution;
    TrilinosWrappers::MPI::Vector     system_rhs;
    parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> soltrans;
#else
    PETScWrappers::MPI::SparseMatrix system_matrix;
    PETScWrappers::MPI::Vector solution;
    PETScWrappers::MPI::Vector system_rhs;
    parallel::distributed::SolutionTransfer<dim, PETScWrappers::MPI::Vector> soltrans;
#endif

    ConditionalOStream                pcout;


};



template <int dim>
LaplaceProblem<dim>::LaplaceProblem ()
		:
		triangulation (MPI_COMM_WORLD
			       ,typename Triangulation<dim>::MeshSmoothing  (Triangulation<dim>::smoothing_on_refinement | Triangulation<dim>::smoothing_on_coarsening)
		),
		dof_handler (triangulation),
                fe (2),
		soltrans(dof_handler),
		pcout (std::cout,
		       (Utilities::System::
			get_this_mpi_process(MPI_COMM_WORLD)
			== 0))
{}



template <int dim>
LaplaceProblem<dim>::~LaplaceProblem ()
{
  dof_handler.clear ();
}



template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  TimeBlock<ConditionalOStream> t(pcout, "setup_system", false);
  {
    TimeBlock<ConditionalOStream> t(pcout, "_init_vectors");
#ifdef USE_TRILINOS
    solution.reinit (locally_relevant_dofs, MPI_COMM_WORLD);
    system_rhs.reinit (locally_owned_dofs, MPI_COMM_WORLD);
#else
    solution.reinit ( MPI_COMM_WORLD, locally_owned_dofs, locally_relevant_dofs);
    system_rhs.reinit (MPI_COMM_WORLD, dof_handler.n_dofs(),
		       dof_handler.n_locally_owned_dofs());

    solution = 0;
    system_rhs = 0;
#endif
  }

  constraints.clear ();
  constraints.reinit (locally_relevant_dofs);

  {

    TimeBlock<ConditionalOStream> t(pcout, "_make_hn");

    DoFTools::make_hanging_node_constraints (static_cast<const DoFHandler<dim>&>(dof_handler),
					     constraints);
  }

  {
    TimeBlock<ConditionalOStream> t(pcout, "_interpol_bv");
    VectorTools::interpolate_boundary_values (static_cast<const DoFHandler<dim>&>(dof_handler),
					      0,
					      ZeroFunction<dim>(),
					      constraints);
    constraints.close ();
  }


  CompressedSimpleSparsityPattern csp (dof_handler.n_dofs(),
				       dof_handler.n_dofs(),
				       locally_relevant_dofs);
  {
    TimeBlock<ConditionalOStream> t(pcout, "_make_sp", false);

    DoFTools::make_sparsity_pattern (static_cast<const DoFHandler<dim>&>(dof_handler),
				     csp,
				     constraints, false);
  }

  {
    TimeBlock<ConditionalOStream> t(pcout, "_init_matrix", false);
#ifdef USE_TRILINOS
    system_matrix.reinit (locally_owned_dofs, csp, MPI_COMM_WORLD, true);
#else

    {
      TimeBlock<ConditionalOStream> t(pcout, "__prep_csp", false);
      SparsityTools::distribute_sparsity_pattern<>(csp,
						   dof_handler.n_locally_owned_dofs_per_processor(),
						   MPI_COMM_WORLD,
						   locally_relevant_dofs);
    }

    system_matrix.reinit (MPI_COMM_WORLD,
			  csp,
			  dof_handler.n_locally_owned_dofs_per_processor(),
			  dof_handler.n_locally_owned_dofs_per_processor(),
			  Utilities::System::get_this_mpi_process(MPI_COMM_WORLD));
#if (PETSC_VERSION_MAJOR <= 2)
//  MatSetOption (system_matrix, MAT_YES_NEW_NONZERO_LOCATIONS);
#else
//  MatSetOption (system_matrix, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE);
//  MatSetOption (system_matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
#endif

#endif
  }
}



template <int dim>
void LaplaceProblem<dim>::assemble_system ()
{
  TimeBlock<ConditionalOStream> t(pcout, "assemble");

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

	      double v = fe_values.quadrature_point(q_point)[0];
	      v = 0.5+0.25*sin(4.0*numbers::PI*v);

	      cell_rhs(i) += (fe_values.shape_value(i,q_point) *
			      (v < fe_values.quadrature_point(q_point)[1]
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

#ifdef USE_TRILINOS
  system_rhs.compress (Add);
#else
  system_rhs.compress ();
#endif
  system_matrix.compress ();
}





template <int dim>
void LaplaceProblem<dim>::solve ()
{
				   // we do not want to get spammed with solver info.
  deallog.depth_console (0);

  TimeBlock<ConditionalOStream> t(pcout, "solve", false);

#ifdef USE_TRILINOS
  TrilinosWrappers::MPI::Vector distributed_solution (locally_owned_dofs,
						      MPI_COMM_WORLD);

  distributed_solution.reinit(solution, false, true);
#else

  PETScWrappers::MPI::Vector distributed_solution(MPI_COMM_WORLD,
						  dof_handler.n_dofs(),
						  dof_handler.n_locally_owned_dofs());

  distributed_solution = solution;
#endif

  print_memory_stats();

  SolverControl solver_control (solution.size(), 1e-12);
#ifdef USE_TRILINOS
  SolverCG<TrilinosWrappers::MPI::Vector> solver (solver_control);

  TrilinosWrappers::PreconditionAMG preconditioner;
  preconditioner.initialize(system_matrix);
#else
  PETScWrappers::SolverCG solver(solver_control, MPI_COMM_WORLD);
  PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);

#endif

  { // memory consumption
    pcout << "Mem: Tria (p4est) DofH Constraints Mat X Rhs I_lo I_lr"
	  << std::endl;
//    MPI_Barrier(MPI_COMM_WORLD);

    pcout << "MEM: "
	  << " " << triangulation.memory_consumption()
	  << " " << triangulation.memory_consumption_p4est()
	  << " " << dof_handler.memory_consumption()
	  << " " << constraints.memory_consumption()
	  << " " << system_matrix.memory_consumption()
	  << " " << solution.memory_consumption()
	  << " " << system_rhs.memory_consumption()
	  << " " << locally_owned_dofs.memory_consumption()
	  << " " << locally_relevant_dofs.memory_consumption()
//	  << " " << preconditioner.memory_consumption()
	  << " ~sum=" << (triangulation.memory_consumption() + dof_handler.memory_consumption()
			  + constraints.memory_consumption() + system_matrix.memory_consumption()
			  + solution.memory_consumption()  + system_rhs.memory_consumption()
			  + locally_owned_dofs.memory_consumption()
			  + locally_relevant_dofs.memory_consumption() )/1024
	  << " kB"
	  << std::endl;

//    MPI_Barrier(MPI_COMM_WORLD);
  }


  solver.solve (system_matrix, distributed_solution, system_rhs,
	    preconditioner);

  print_memory_stats();

  pcout << "  Solved in " << solver_control.last_step()
	<< " iterations." << std::endl;

  constraints.distribute (distributed_solution);
  solution = distributed_solution;

#ifndef USE_TRILINOS
  solution.update_ghost_values();
#endif

  deallog.depth_console (1);
}



template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
  TimerOutput computing_timer(pcout, TimerOutput::summary,
			      TimerOutput::wall_times);

  computing_timer.enter_section ("Estimating");

  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  {

    TimeBlock<ConditionalOStream> t(pcout, "kelly");

  KellyErrorEstimator<dim>::estimate (static_cast<const DoFHandler<dim>&>(dof_handler),
				      QGauss<dim-1>(3),
				      typename FunctionMap<dim>::type(),
				      solution,
 				      estimated_error_per_cell);
  }

  computing_timer.exit_section();

  computing_timer.enter_section ("Marking");
  {
    TimeBlock<ConditionalOStream> t(pcout, "marking");

  parallel::distributed::GridRefinement::
    refine_and_coarsen_fixed_number (triangulation,
				     estimated_error_per_cell,
				     0.3, 0.03);
  }

  computing_timer.exit_section();

  computing_timer.enter_section ("Refining");


  {
    TimeBlock<ConditionalOStream> t(pcout, "prep_c&r");

    triangulation.prepare_coarsening_and_refinement();
  }

  if (0)
    {
				       //output refinement information for
				       //debugging

     const std::string filename = ("ref." +
				   Utilities::int_to_string
				   (triangulation.locally_owned_subdomain(), 4) +
				   ".vtk");
     std::ofstream output (filename.c_str());

     DataOut<dim> data_out;
     data_out.attach_dof_handler (dof_handler);

     unsigned int n_coarse=0;
     Vector<float> subdomain (triangulation.n_active_cells());
     {
       unsigned int index = 0;

       for (typename Triangulation<dim>::active_cell_iterator
	      cell = triangulation.begin_active();
	    cell != triangulation.end(); ++cell, ++index)
	 {
	   subdomain(index)=0;

	   if (cell->is_ghost() || cell->is_artificial())
	     subdomain(index)=-4;

	   if (cell->refine_flag_set())
	     subdomain(index)+=1;
	   if (cell->coarsen_flag_set())
	     {
	       subdomain(index)+=2;
	       ++n_coarse;
	     }
	 }
     }
     std::cout << "id=" << triangulation.locally_owned_subdomain() << " n_coarsen=" << n_coarse << std::endl;

     data_out.add_data_vector (subdomain, "info");
     data_out.build_patches ();
     data_out.write_vtk (output);
    }

  {
    TimeBlock<ConditionalOStream> t(pcout, "soltrans.prep");

  soltrans.prepare_for_coarsening_and_refinement(solution);
  }

  {
    TimeBlock<ConditionalOStream> t(pcout, "c&r", false);

  triangulation.execute_coarsening_and_refinement ();
  }

  computing_timer.exit_section();
}




template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
  const std::string filename = ("solution-" +
				Utilities::int_to_string (cycle, 2) +
				"." +
				Utilities::int_to_string
				(triangulation.locally_owned_subdomain(), 4) +
				".d2");
  std::ofstream output (filename.c_str());

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "u");

  Vector<float> subdomain (triangulation.n_active_cells());
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

  data_out.write_deal_II_intermediate (output);
}



template <int dim>
void LaplaceProblem<dim>::run (unsigned int refine)
{
      {
	TimeBlock<ConditionalOStream> t(pcout, "empty");
      }

      unsigned int cycles = 12;

      pcout << "running " <<  cycles << " cycles with refine= " << refine <<
#ifdef USE_TRILINOS
	    " with Trilinos"
#else
	    " with PETSc"
#endif
	    << std::endl;

  for (unsigned int cycle=0; cycle<cycles; ++cycle)
    {
      TimerOutput computing_timer(pcout, TimerOutput::summary,
				  TimerOutput::wall_times);

      pcout << "Cycle " << cycle << ':' << std::endl;

      computing_timer.enter_section ("Mesh handling");
      if (cycle == 0)
	{
	  GridGenerator::hyper_cube (triangulation);
	  triangulation.refine_global (refine);
	}
      else
	refine_grid ();



      pcout << "   Number of active cells:       "
	    << triangulation.n_global_active_cells()
	    << std::endl;
/*      pcout << "                                 (";
      for (unsigned int i=0; i<triangulation.n_active_cells().size(); ++i)
	pcout << triangulation.n_active_cells()[i]
	      << (i != triangulation.n_active_cells().size()-1 ? "+" : "");
	      pcout << ")" << std::endl;*/
      computing_timer.exit_section();

      computing_timer.enter_section ("DoFs handling");

      {
	TimeBlock<ConditionalOStream> t(pcout, "dist_dofs");
	dof_handler.distribute_dofs (fe);
      }


      locally_owned_dofs = dof_handler.locally_owned_dofs ();
      DoFTools::extract_locally_relevant_dofs (dof_handler,
					       locally_relevant_dofs);

      pcout << "   Number of degrees of freedom: "
	    << dof_handler.n_dofs()
	    << std::endl;
/*      pcout << "                                 (";
      for (unsigned int i=0; i<dof_handler.n_dofs().size(); ++i)
	pcout << dof_handler.n_dofs()[i]
	      << (i != dof_handler.n_dofs().size()-1 ? "+" : "");
	      pcout << ")" << std::endl;*/
      computing_timer.exit_section();

      computing_timer.enter_section ("Setting up system");
      setup_system ();

      if (cycle>0)
	{
	  TimeBlock<ConditionalOStream> t(pcout, "soltrans.interp");

	soltrans.interpolate(solution);
	}


      computing_timer.exit_section();

      computing_timer.enter_section ("Assembling system");
      assemble_system ();
      computing_timer.exit_section();

      computing_timer.enter_section ("Solving system");
      solve ();
      computing_timer.exit_section();

//      output_results (cycle);

      pcout << std::endl;
    }
}



int main(int argc, char *argv[])
{
  try
    {
#ifdef USE_TRILINOS
      Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv);
#else
      PetscInitialize(&argc, &argv);
#endif

      deallog.depth_console (1);

      int refine=5;
      if (argc>1)
	{
	  refine = (unsigned int)Utilities::string_to_int(argv[1]);
	}

      {
	LaplaceProblem<2> laplace_problem_2d;
	laplace_problem_2d.run (refine);
      }

      print_memory_stats();

#ifndef USE_TRILINOS
      PetscLogPrintSummary(MPI_COMM_WORLD, "petsc.log" );
      PetscFinalize();
#endif

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
