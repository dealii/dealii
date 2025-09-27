// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// This test provides a minimal Laplace equation in 2D for testing
// the BDDC preconditioner. It cannot be simpler as it fails to setup
// in simpler matrices (i.e. diagonal), and the preset matrix factories
// create non-IS matrices, which are not supported by BDDC.

#include <deal.II/base/index_set.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/linear_operator.h>

#include "../tests.h"

// Vectors:
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

// Block Matrix and Vectors:
#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_block_vector.h>

// Dof and sparsity tools:
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_tools.h>



int
main(int argc, char *argv[])
{
  using size_type = PETScWrappers::MPI::SparseMatrix::size_type;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  deallog << std::setprecision(10);
  parallel::distributed::Triangulation<2> triangulation(MPI_COMM_WORLD);
  FE_Q<2>                                 fe(1);
  DoFHandler<2>                           dof_handler(triangulation);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  dof_handler.distribute_dofs(fe);

  const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet  locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);
  const IndexSet locally_active_dofs =
    DoFTools::extract_locally_active_dofs(dof_handler);

  DynamicSparsityPattern dsp(locally_relevant_dofs);

  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  SparsityTools::distribute_sparsity_pattern(dsp,
                                             dof_handler.locally_owned_dofs(),
                                             MPI_COMM_WORLD,
                                             locally_relevant_dofs);

  PETScWrappers::MPI::SparseMatrix system_matrix;
  system_matrix.reinit(locally_owned_dofs,
                       locally_active_dofs,
                       locally_owned_dofs,
                       locally_active_dofs,
                       dsp,
                       MPI_COMM_WORLD);
  deallog << "MATIS:OK" << std::endl;
  PETScWrappers::MPI::Vector locally_relevant_solution(locally_owned_dofs,
                                                       locally_relevant_dofs,
                                                       MPI_COMM_WORLD);
  PETScWrappers::MPI::Vector completely_distributed_solution(locally_owned_dofs,
                                                             MPI_COMM_WORLD);
  PETScWrappers::MPI::Vector system_rhs(locally_owned_dofs, MPI_COMM_WORLD);

  const QGauss<2> quadrature_formula(2);
  FEValues<2>     fe_values(fe,
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
        cell_matrix = 0.;
        cell_rhs    = 0.;

        fe_values.reinit(cell);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const double rhs_value = 1.0;

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  {
                    cell_matrix(i, j) += fe_values.shape_grad(i, q_point) *
                                         fe_values.shape_grad(j, q_point) *
                                         fe_values.JxW(q_point);
                    cell_matrix(i, j) += fe_values.shape_value(i, q_point) *
                                         fe_values.shape_value(j, q_point) *
                                         fe_values.JxW(q_point);
                  }

                cell_rhs(i) += 1.0 * fe_values.shape_value(i, q_point) *
                               fe_values.JxW(q_point);
              }
          }

        cell->get_dof_indices(local_dof_indices);
        system_matrix.add(local_dof_indices, cell_matrix);
        system_rhs.add(local_dof_indices, cell_rhs);
      }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

  SolverControl solver_control(dof_handler.n_dofs(), 1e-12);

  PETScWrappers::SolverCG                            solver(solver_control);
  PETScWrappers::PreconditionBDDC<2>                 preconditioner;
  PETScWrappers::PreconditionBDDC<2>::AdditionalData data;

  // Now we setup the dof coordinates if a sufficiently new PETSc is used
  std::map<types::global_dof_index, Point<2>> dof_2_point;
  DoFTools::map_dofs_to_support_points(MappingQ1<2>(),
                                       dof_handler,
                                       dof_2_point);
  std::vector<Point<2>> coords(locally_owned_dofs.n_elements());
  unsigned int          k = 0;
  for (auto it = locally_owned_dofs.begin(); it != locally_owned_dofs.end();
       ++it, ++k)
    {
      coords[k] = dof_2_point[*it];
    }
  data.coords       = coords;
  data.use_vertices = true;

  preconditioner.initialize(system_matrix, data);
  check_solver_within_range(solver.solve(system_matrix,
                                         completely_distributed_solution,
                                         system_rhs,
                                         preconditioner),
                            solver_control.last_step(),
                            1,
                            2);

  deallog << "CG/BDDC:OK" << std::endl;
  return 0;
}
