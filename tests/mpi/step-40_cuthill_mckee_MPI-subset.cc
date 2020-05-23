// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// This is a variation of the step-40_cuthill_mckee test, but instead
// of working on MPI_COMM_WORLD, it works on only a subset of
// processes. The idea of the test is to verify that we can do that,
// i.e., that we don't have any accidental references to
// MPI_COMM_WORLD that would lead to deadlocks if we were working on a
// different communicator.

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

namespace Step40
{
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem(MPI_Comm comm);
    ~LaplaceProblem();

    void
    run();

  private:
    void
    setup_system();
    void
    assemble_system();
    void
    solve();
    void
    refine_grid();

    MPI_Comm mpi_communicator;

    parallel::distributed::Triangulation<dim> triangulation;

    DoFHandler<dim> dof_handler;
    FE_Q<dim>       fe;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    AffineConstraints<double> constraints;

    PETScWrappers::MPI::SparseMatrix system_matrix;
    PETScWrappers::MPI::Vector       locally_relevant_solution;
    PETScWrappers::MPI::Vector       system_rhs;

    ConditionalOStream pcout;
  };



  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem(MPI_Comm comm)
    : mpi_communicator(comm)
    , triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening))
    , dof_handler(triangulation)
    , fe(2)
    , pcout(Utilities::MPI::this_mpi_process(mpi_communicator) == 0 ?
              deallog.get_file_stream() :
              std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  {}



  template <int dim>
  LaplaceProblem<dim>::~LaplaceProblem()
  {
    dof_handler.clear();
  }



  template <int dim>
  void
  LaplaceProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    {
      std::vector<types::global_dof_index> starting_indices;


      const QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);
      FEFaceValues<dim>     fe_face_values(fe,
                                       face_quadrature_formula,
                                       update_values |
                                         update_quadrature_points |
                                         update_normal_vectors |
                                         update_JxW_values);

      Tensor<1, dim>                       u;
      Point<dim>                           down{0, -1};
      std::vector<types::global_dof_index> dof_indices(fe.n_dofs_per_face(), 0);


      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              for (const unsigned int face : GeometryInfo<dim>::face_indices())
                {
                  if ((cell->face(face)->at_boundary()) ||
                      (cell->neighbor(face)->is_active() &&
                       cell->neighbor(face)->is_ghost()))
                    {
                      fe_face_values.reinit(cell, face);
                      // for Q_2 this is in middle of face, dim=2 or what
                      // quadrature point to give?
                      u = fe_face_values.normal_vector(1);
                      if (u * down < 0)
                        {
                          cell->face(face)->get_dof_indices(dof_indices);
                          starting_indices.insert(std::end(starting_indices),
                                                  std::begin(dof_indices),
                                                  std::end(dof_indices));
                        }
                    }
                }
            }
        }

      // remove duplicates by creating a set
      std::set<types::global_dof_index> no_duplicates_please(
        starting_indices.begin(), starting_indices.end());
      // back to vector for the DoFRenumbering function
      starting_indices.clear();
      starting_indices.assign(no_duplicates_please.begin(),
                              no_duplicates_please.end());

      DoFRenumbering::Cuthill_McKee(dof_handler, false, true, starting_indices);

      locally_owned_dofs = dof_handler.locally_owned_dofs();
      DoFTools::extract_locally_relevant_dofs(dof_handler,
                                              locally_relevant_dofs);
    }

    locally_relevant_solution.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    system_rhs = PetscScalar();

    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             constraints);
    constraints.close();

    DynamicSparsityPattern csp(dof_handler.n_dofs(),
                               dof_handler.n_dofs(),
                               locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, csp, constraints, false);
    SparsityTools::distribute_sparsity_pattern(csp,
                                               locally_owned_dofs,
                                               mpi_communicator,
                                               locally_relevant_dofs);
    system_matrix.reinit(
      mpi_communicator,
      csp,
      Utilities::MPI::all_gather(mpi_communicator,
                                 dof_handler.n_locally_owned_dofs()),
      Utilities::MPI::all_gather(mpi_communicator,
                                 dof_handler.n_locally_owned_dofs()),
      Utilities::MPI::this_mpi_process(mpi_communicator));
  }



  template <int dim>
  void
  LaplaceProblem<dim>::assemble_system()
  {
    const QGauss<dim> quadrature_formula(3);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<PetscScalar> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<PetscScalar>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
        {
          cell_matrix = PetscScalar();
          cell_rhs    = PetscScalar();

          fe_values.reinit(cell);

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
              const double rhs_value =
                (fe_values.quadrature_point(q_point)[1] >
                     0.5 +
                       0.25 * std::sin(4.0 * numbers::PI *
                                       fe_values.quadrature_point(q_point)[0]) ?
                   1 :
                   -1);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    cell_matrix(i, j) += (fe_values.shape_grad(i, q_point) *
                                          fe_values.shape_grad(j, q_point) *
                                          fe_values.JxW(q_point));

                  cell_rhs(i) +=
                    (rhs_value * fe_values.shape_value(i, q_point) *
                     fe_values.JxW(q_point));
                }
            }

          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 cell_rhs,
                                                 local_dof_indices,
                                                 system_matrix,
                                                 system_rhs);
        }

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }



  template <int dim>
  void
  LaplaceProblem<dim>::solve()
  {
    PETScWrappers::MPI::Vector completely_distributed_solution(
      mpi_communicator,
      dof_handler.n_dofs(),
      dof_handler.n_locally_owned_dofs());

    SolverControl solver_control(dof_handler.n_dofs(), 1e-12);

    PETScWrappers::SolverCG solver(solver_control, mpi_communicator);

#ifndef PETSC_USE_COMPLEX
    PETScWrappers::PreconditionBoomerAMG preconditioner(
      system_matrix,
      PETScWrappers::PreconditionBoomerAMG::AdditionalData(true));

    check_solver_within_range(solver.solve(system_matrix,
                                           completely_distributed_solution,
                                           system_rhs,
                                           preconditioner),
                              solver_control.last_step(),
                              11,
                              11);
#else
    check_solver_within_range(solver.solve(system_matrix,
                                           completely_distributed_solution,
                                           system_rhs,
                                           PETScWrappers::PreconditionJacobi(
                                             system_matrix)),
                              solver_control.last_step(),
                              120,
                              260);
#endif

    pcout << "   Solved in " << solver_control.last_step() << " iterations."
          << std::endl;

    constraints.distribute(completely_distributed_solution);

    locally_relevant_solution = completely_distributed_solution;
  }



  template <int dim>
  void
  LaplaceProblem<dim>::refine_grid()
  {
    triangulation.refine_global(1);
  }



  template <int dim>
  void
  LaplaceProblem<dim>::run()
  {
    const unsigned int n_cycles = 2;
    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
      {
        pcout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube(triangulation);
            triangulation.refine_global(5);
          }
        else
          refine_grid();

        setup_system();

        pcout << "   Number of active cells:       "
              << triangulation.n_global_active_cells() << std::endl
              << "      ";
        const auto n_locally_owned_active_cells_per_processor =
          Utilities::MPI::all_gather(
            triangulation.get_communicator(),
            triangulation.n_locally_owned_active_cells());
        for (unsigned int i = 0;
             i < Utilities::MPI::n_mpi_processes(mpi_communicator);
             ++i)
          pcout << n_locally_owned_active_cells_per_processor[i] << '+';
        pcout << std::endl;

        pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl
              << "      ";
        const std::vector<types::global_dof_index>
          n_locally_owned_dofs_per_processor =
            Utilities::MPI::all_gather(mpi_communicator,
                                       dof_handler.n_locally_owned_dofs());
        for (unsigned int i = 0;
             i < Utilities::MPI::n_mpi_processes(mpi_communicator);
             ++i)
          pcout << n_locally_owned_dofs_per_processor[i] << '+';
        pcout << std::endl;

        assemble_system();
        solve();

        pcout << std::endl;
      }
  }
} // namespace Step40


int
test_mpi(MPI_Comm comm)
{
  try
    {
      using namespace Step40;


      {
        LaplaceProblem<2> laplace_problem_2d(comm);
        laplace_problem_2d.run();
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
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
      std::cerr << std::endl
                << std::endl
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



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  // create a group of 4 out of the 7 processes that we want to run
  // this program with
  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  Assert(n_procs == 7, ExcInternalError());

  MPI_Group whole_group;
  MPI_Comm_group(MPI_COMM_WORLD, &whole_group);

  MPI_Group subset_group;
  const int subset_ranks[] = {6, 0, 2, 3};
  MPI_Group_incl(whole_group, 4, subset_ranks, &subset_group);

  MPI_Comm subset_comm;
  MPI_Comm_create(MPI_COMM_WORLD, subset_group, &subset_comm);

  // now only run the program on the subset of processors identified
  // in 'comm'. all of the other processes simply do nothing and will
  // wait for termination in the destructor of MPI_InitFinalize until
  // the worker processes are ready to join them
  if (std::find(std::begin(subset_ranks),
                std::end(subset_ranks),
                Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)) !=
      std::end(subset_ranks))
    {
      if (Utilities::MPI::this_mpi_process(subset_comm) == 0)
        {
          initlog();

          // check that creation above worked correctly
          Assert(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 6,
                 ExcInternalError());
          Assert(Utilities::MPI::n_mpi_processes(subset_comm) == 4,
                 ExcInternalError());

          deallog.push("mpi");
          test_mpi(subset_comm);
          deallog.pop();
        }
      else
        test_mpi(subset_comm);
    }
}
