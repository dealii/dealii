// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

// test VectorTools::interpolate_to_different_mesh in parallel
// this is a slightly modified version from the example by Sam Cox from the
// mailing


namespace LA
{
  using namespace dealii::LinearAlgebraTrilinos;
}

template <int dim>
class SeventhProblem
{
public:
  SeventhProblem(unsigned int prob_number);
  ~SeventhProblem();
  void
  run(unsigned int cycle);

private:
  void
  setup_system();
  void
  setup_second_system();
  void
  assemble_system();
  void
                                                    solve();
  MPI_Comm                                          mpi_communicator;
  parallel::distributed::Triangulation<2>::Settings settings;
  parallel::distributed::Triangulation<dim>         triangulation;
  DoFHandler<dim>                                   dof_handler;
  FE_Q<dim>                                         fe;
  IndexSet                                          locally_owned_dofs;
  IndexSet                                          locally_relevant_dofs;
  AffineConstraints<double>                         constraints;
  TrilinosWrappers::SparseMatrix                    system_matrix;
  TrilinosWrappers::MPI::Vector                     locally_relevant_solution;
  TrilinosWrappers::MPI::Vector interpolated_locally_relevant_solution;
  TrilinosWrappers::MPI::Vector system_rhs;
  parallel::distributed::Triangulation<dim> second_triangulation;
  DoFHandler<dim>                           second_dof_handler;
  FE_Q<dim>                                 second_fe;
  IndexSet                                  second_locally_owned_dofs;
  IndexSet                                  second_locally_relevant_dofs;
  TrilinosWrappers::MPI::Vector             second_locally_relevant_solution;
  ConditionalOStream                        pcout;
  unsigned int                              prob_number;
};
template <int dim>
SeventhProblem<dim>::SeventhProblem(unsigned int prob_number)
  : mpi_communicator(MPI_COMM_WORLD)
  , settings(
      parallel::distributed::Triangulation<2, 2>::no_automatic_repartitioning)
  , triangulation(mpi_communicator,
                  typename Triangulation<dim>::MeshSmoothing(
                    Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening),
                  settings)
  , dof_handler(triangulation)
  , fe(2)
  , second_triangulation(mpi_communicator,
                         typename Triangulation<dim>::MeshSmoothing(
                           Triangulation<dim>::smoothing_on_refinement |
                           Triangulation<dim>::smoothing_on_coarsening),
                         settings)
  , second_dof_handler(second_triangulation)
  , second_fe(2)
  , pcout(deallog.get_file_stream(),
          (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  , prob_number(prob_number)
{}

template <int dim>
SeventhProblem<dim>::~SeventhProblem()
{
  dof_handler.clear();
  second_dof_handler.clear();
}

template <int dim>
void
SeventhProblem<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);
  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
  locally_relevant_solution.reinit(locally_owned_dofs,
                                   locally_relevant_dofs,
                                   mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
  system_rhs = 0;
  constraints.clear();
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();
  DynamicSparsityPattern csp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, csp, constraints, false);
  SparsityTools::distribute_sparsity_pattern(csp,
                                             locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       csp,
                       mpi_communicator);
}

template <int dim>
void
SeventhProblem<dim>::setup_second_system()
{
  second_dof_handler.distribute_dofs(fe);
  second_locally_owned_dofs = second_dof_handler.locally_owned_dofs();
  second_locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(second_dof_handler);
  second_locally_relevant_solution.reinit(second_locally_owned_dofs,
                                          second_locally_relevant_dofs,
                                          mpi_communicator);
  interpolated_locally_relevant_solution.reinit(second_locally_owned_dofs,
                                                second_locally_relevant_dofs,
                                                mpi_communicator);
}

template <int dim>
void
SeventhProblem<dim>::assemble_system()
{
  const QGauss<dim>  quadrature_formula(3);
  FEValues<dim>      fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    if (cell->is_locally_owned())
      {
        cell_matrix = 0;
        cell_rhs    = 0;
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
                  cell_matrix(i, j) +=
                    (fe_values.shape_grad(i, q_point) *
                     fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point));
                cell_rhs(i) += (rhs_value * fe_values.shape_value(i, q_point) *
                                fe_values.JxW(q_point));
              }
          }
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}
template <int dim>
void
SeventhProblem<dim>::solve()
{
  LA::MPI::Vector            completely_distributed_solution(locally_owned_dofs,
                                                  mpi_communicator);
  SolverControl              solver_control(dof_handler.n_dofs(), 1e-12);
  TrilinosWrappers::SolverCG solver(solver_control);
  TrilinosWrappers::PreconditionAMG                 preconditioner;
  TrilinosWrappers::PreconditionAMG::AdditionalData data;
  preconditioner.initialize(system_matrix, data);
  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
               preconditioner);
  pcout << " Solved in " << solver_control.last_step() << " iterations."
        << std::endl;
  constraints.distribute(completely_distributed_solution);
  locally_relevant_solution = completely_distributed_solution;
}

template <int dim>
void
SeventhProblem<dim>::run(unsigned int cycle)
{
  if (cycle == 0)
    {
      GridGenerator::hyper_cube(triangulation);
      triangulation.refine_global(1);
      GridGenerator::hyper_cube(second_triangulation);
      second_triangulation.refine_global(1);
      setup_system();
    }
  else
    {
      Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
      KellyErrorEstimator<dim>::estimate(
        dof_handler,
        QGauss<dim - 1>(3),
        std::map<types::boundary_id, const Function<dim> *>(),
        locally_relevant_solution,
        estimated_error_per_cell);

      parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
        triangulation, estimated_error_per_cell, 0.5, 0.3);
      std::vector<bool> r_flags;
      std::vector<bool> c_flags;
      triangulation.prepare_coarsening_and_refinement();
      triangulation.save_refine_flags(r_flags);
      triangulation.save_coarsen_flags(c_flags);

      triangulation.execute_coarsening_and_refinement();

      setup_system();
      pcout << " Number of active cells: "
            << triangulation.n_global_active_cells() << std::endl
            << " Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
      assemble_system();
      solve();

      setup_second_system();
      second_locally_relevant_solution = locally_relevant_solution;

      VectorTools::interpolate_to_different_mesh(
        dof_handler,
        locally_relevant_solution,
        second_dof_handler,
        interpolated_locally_relevant_solution);
      second_triangulation.load_coarsen_flags(c_flags);
      second_triangulation.load_refine_flags(r_flags);
      second_triangulation.execute_coarsening_and_refinement();
    }
}

void
seventh_grid()
{
  ConditionalOStream pcout(deallog.get_file_stream(),
                           (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                            0));

  pcout << "7th Starting" << std::endl;
  SeventhProblem<2>  lap(1);
  const unsigned int n_cycles = 5;
  for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
    {
      pcout << "Cycle " << cycle << ':' << std::endl;
      lap.run(cycle);
    }
  pcout << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  deallog.depth_file(0);

  seventh_grid();
}
