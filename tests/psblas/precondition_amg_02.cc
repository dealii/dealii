// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2026 by the deal.II authors
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
#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/types.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/psblas_precondition.h>
#include <deal.II/lac/psblas_sparse_matrix.h>
#include <deal.II/lac/psblas_vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

using namespace dealii;

// Solve a Laplace problem with CG preconditioned by AMG4PSBLAS on a sequence
// of adaptively refined meshes, as in step-40.

template <int dim>
class LaplaceProblem
{
public:
  LaplaceProblem();

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

  const FE_Q<dim> fe;
  DoFHandler<dim> dof_handler;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  AffineConstraints<double> constraints;

  PSCToolkitWrappers::SparseMatrix system_matrix;
  PSCToolkitWrappers::Vector       locally_relevant_solution;
  PSCToolkitWrappers::Vector       system_rhs;
};



template <int dim>
LaplaceProblem<dim>::LaplaceProblem()
  : mpi_communicator(MPI_COMM_WORLD)
  , triangulation(mpi_communicator,
                  typename Triangulation<dim>::MeshSmoothing(
                    Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening))
  , fe(1)
  , dof_handler(triangulation)
{}



template <int dim>
void
LaplaceProblem<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  locally_relevant_solution.reinit(locally_owned_dofs,
                                   locally_relevant_dofs,
                                   mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);

  constraints.clear();
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           types::boundary_id(0),
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  system_matrix.reinit(locally_owned_dofs, mpi_communicator);
}



template <int dim>
void
LaplaceProblem<dim>::assemble_system()
{
  const QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);

        cell_matrix = 0.;
        cell_rhs    = 0.;

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const double rhs_value =
              (fe_values.quadrature_point(q_point)[1] >
                   0.5 +
                     0.25 * std::sin(4.0 * numbers::PI *
                                     fe_values.quadrature_point(q_point)[0]) ?
                 1. :
                 -1.);

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) += fe_values.shape_grad(i, q_point) *
                                       fe_values.shape_grad(j, q_point) *
                                       fe_values.JxW(q_point);

                cell_rhs(i) += rhs_value * fe_values.shape_value(i, q_point) *
                               fe_values.JxW(q_point);
              }
          }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }

  system_matrix.compress();
  system_rhs.compress(VectorOperation::add);
}



template <int dim>
void
LaplaceProblem<dim>::solve()
{
  PSCToolkitWrappers::Vector completely_distributed_solution(locally_owned_dofs,
                                                             mpi_communicator);

  SolverControl                        solver_control(dof_handler.n_dofs(),
                               1e-8 * system_rhs.l2_norm(),
                               false,
                               false);
  SolverCG<PSCToolkitWrappers::Vector> solver(solver_control);

  PSCToolkitWrappers::PreconditionAMG                 preconditioner;
  PSCToolkitWrappers::PreconditionAMG::AdditionalData prec_data;
  prec_data.cycle_type  = "VCYCLE";
  prec_data.aggr_prol   = "SMOOTHED";
  prec_data.n_cycles    = 1;
  prec_data.coarse_type = "ILU";
  preconditioner.initialize(system_matrix, prec_data);

  check_solver_within_range(solver.solve(system_matrix,
                                         completely_distributed_solution,
                                         system_rhs,
                                         preconditioner),
                            solver_control.last_step(),
                            1,
                            40);

  const double norm_before = completely_distributed_solution.l2_norm();
  AssertThrow(norm_before > 0.,
              ExcMessage("The solver returned a zero solution."));


  constraints.distribute(completely_distributed_solution);

  const double norm_after = completely_distributed_solution.l2_norm();
  AssertThrow(norm_after > 0.5 * norm_before,
              ExcMessage("AffineConstraints::distribute() discarded parts of "
                         "the solution."));

  locally_relevant_solution = completely_distributed_solution;

  AssertThrow(locally_relevant_solution.l2_norm() > 0.,
              ExcMessage("The ghosted solution vector is zero."));
}



template <int dim>
void
LaplaceProblem<dim>::refine_grid()
{
  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
  KellyErrorEstimator<dim>::estimate(
    dof_handler,
    QGauss<dim - 1>(fe.degree + 1),
    std::map<types::boundary_id, const Function<dim> *>(),
    locally_relevant_solution,
    estimated_error_per_cell);

  AssertThrow(estimated_error_per_cell.linfty_norm() > 0.,
              ExcMessage("The error estimator returned zero everywhere, which "
                         "means that the solution was lost."));

  // Coarsening is switched off so that the number of cells can only grow.
  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
    triangulation, estimated_error_per_cell, 0.3, 0.0);
  triangulation.execute_coarsening_and_refinement();
}



template <int dim>
void
LaplaceProblem<dim>::run()
{
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  types::global_cell_index previous_n_cells = 0;

  for (unsigned int cycle = 0; cycle < 4; ++cycle)
    {
      deallog << "Cycle " << cycle << std::endl;

      if (cycle > 0)
        refine_grid();

      AssertThrow(triangulation.n_global_active_cells() > previous_n_cells,
                  ExcMessage("The mesh did not grow during refinement."));
      previous_n_cells = triangulation.n_global_active_cells();

      setup_system();
      assemble_system();
      solve();
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 2,
              ExcMessage("This test needs to be run with 2 MPI processes."));

  initlog();

  LaplaceProblem<2> laplace_problem;
  laplace_problem.run();

  deallog << "OK" << std::endl;

  return 0;
}
