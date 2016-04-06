// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// A lightly adapted version of the step-40 tutorial program

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <fstream>
#include <iostream>

namespace Step40
{
  using namespace dealii;


  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem ();
    ~LaplaceProblem ();

    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();

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
    pcout (Utilities::MPI::this_mpi_process(mpi_communicator)
           == 0
           ?
           deallog.get_file_stream()
           :
           std::cout,
           (Utilities::MPI::this_mpi_process(mpi_communicator)
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
    dof_handler.distribute_dofs (fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler,
                                             locally_relevant_dofs);

    locally_relevant_solution.reinit (locally_owned_dofs,
                                      locally_relevant_dofs,
                                      mpi_communicator);
    system_rhs.reinit (locally_owned_dofs, mpi_communicator);
    system_rhs = PetscScalar();

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
                          Utilities::MPI::this_mpi_process(mpi_communicator));
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

    FullMatrix<PetscScalar>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<PetscScalar>       cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          cell_matrix = PetscScalar();
          cell_rhs = PetscScalar();

          fe_values.reinit (cell);

          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            {
              const double
              rhs_value
                = (fe_values.quadrature_point(q_point)[1]
                   >
                   0.5+0.25*std::sin(4.0 * numbers::PI *
                                     fe_values.quadrature_point(q_point)[0])
                   ? 1 : -1);

              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
                                         fe_values.shape_grad(j,q_point) *
                                         fe_values.JxW(q_point));

                  cell_rhs(i) += (rhs_value *
                                  fe_values.shape_value(i,q_point) *
                                  fe_values.JxW(q_point));
                }
            }

          cell->get_dof_indices (local_dof_indices);
          constraints.distribute_local_to_global (cell_matrix,
                                                  cell_rhs,
                                                  local_dof_indices,
                                                  system_matrix,
                                                  system_rhs);
        }

    system_matrix.compress (VectorOperation::add);
    system_rhs.compress (VectorOperation::add);
  }




  template <int dim>
  void LaplaceProblem<dim>::solve ()
  {
    PETScWrappers::MPI::Vector
    completely_distributed_solution (mpi_communicator,
                                     dof_handler.n_dofs(),
                                     dof_handler.n_locally_owned_dofs());

    SolverControl solver_control (dof_handler.n_dofs(), 1e-12);
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
    solver.set_symmetric_mode(true);
    solver.solve(system_matrix, completely_distributed_solution, system_rhs);

    pcout << "   Solved in " << solver_control.last_step()
          << " iterations." << std::endl;

    constraints.distribute (completely_distributed_solution);

    locally_relevant_solution = completely_distributed_solution;
  }




  template <int dim>
  void LaplaceProblem<dim>::refine_grid ()
  {
    triangulation.refine_global (1);
  }




  template <int dim>
  void LaplaceProblem<dim>::run ()
  {
    const unsigned int n_cycles = 2;
    for (unsigned int cycle=0; cycle<n_cycles; ++cycle)
      {
        pcout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube (triangulation);
            triangulation.refine_global (5);
          }
        else
          refine_grid ();

        setup_system ();

        pcout << "   Number of active cells:       "
              << triangulation.n_global_active_cells()
              << std::endl
              << "      ";
        for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(mpi_communicator); ++i)
          pcout << triangulation.n_locally_owned_active_cells_per_processor()[i]
                << '+';
        pcout << std::endl;

        pcout << "   Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl
              << "      ";
        for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(mpi_communicator); ++i)
          pcout << dof_handler.n_locally_owned_dofs_per_processor()[i]
                << '+';
        pcout << std::endl;

        assemble_system ();
        solve ();

        pcout << std::endl;
      }
  }
}


int test_mpi ()
{

  try
    {
      using namespace dealii;
      using namespace Step40;


      {
        LaplaceProblem<2> laplace_problem_2d;
        laplace_problem_2d.run ();
      }
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




int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.threshold_double(1.e-10);

      deallog.push("mpi");
      test_mpi();
      deallog.pop();
    }
  else
    test_mpi();
}
