// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check different functions from the SparseMatrixTools namespace

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix_tools.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/numerics/matrix_creator.h>

#include "../tests.h"

template <int dim, int spacedim>
void
reinit_sparsity_pattern(const DoFHandler<dim, spacedim> &dof_handler,
                        SparsityPattern                 &sparsity_pattern)
{
  std::vector<unsigned int> counter(dof_handler.n_dofs(), 0);

  std::vector<types::global_dof_index> local_dof_indices;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      local_dof_indices.resize(cell->get_fe().n_dofs_per_cell());

      for (const auto i : local_dof_indices)
        counter[i]++;
    }

  sparsity_pattern.reinit(dof_handler.n_dofs(),
                          dof_handler.n_dofs(),
                          *std::max_element(counter.begin(), counter.end()));
}

template <int dim, int spacedim>
void
reinit_sparsity_pattern(const DoFHandler<dim, spacedim>   &dof_handler,
                        TrilinosWrappers::SparsityPattern &sparsity_pattern)
{
  sparsity_pattern.reinit(dof_handler.locally_owned_dofs(),
                          dof_handler.get_mpi_communicator());
}

template <int dim,
          typename SpareMatrixType,
          typename SparsityPatternType,
          typename SparseMatrixType2,
          typename SparsityPatternType2>
void
test()
{
  const unsigned int fe_degree = 1;

  // create mesh, ...
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_cube(tria, 3);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FE_Q<dim>(fe_degree));

  QGauss<dim> quadrature(fe_degree + 1);

  AffineConstraints<double> constraints;
  DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
  constraints.close();

  // create system matrix
  SparsityPatternType sparsity_pattern;
  reinit_sparsity_pattern(dof_handler, sparsity_pattern);
  DoFTools::make_sparsity_pattern(dof_handler,
                                  sparsity_pattern,
                                  constraints,
                                  false);
  sparsity_pattern.compress();

  SpareMatrixType laplace_matrix;
  laplace_matrix.reinit(sparsity_pattern);

  MatrixCreator::create_laplace_matrix<dim, dim, SpareMatrixType>(
    dof_handler, quadrature, laplace_matrix, nullptr, constraints);

  // extract blocks
  std::vector<FullMatrix<double>> blocks;
  SparseMatrixTools::restrict_to_cells(laplace_matrix,
                                       sparsity_pattern,
                                       dof_handler,
                                       blocks);

  for (const auto &block : blocks)
    {
      if (block.m() == 0 && block.m() == 0)
        continue;

      block.print_formatted(deallog.get_file_stream(), 2, false, 8);
      deallog << std::endl;
    }

  const auto test_restrict = [&](const IndexSet &is_0, const IndexSet &is_1) {
    (void)is_1;
    SparsityPatternType2 serial_sparsity_pattern;
    SparseMatrixType2    serial_sparse_matrix;

    if (is_1.size() == 0)
      SparseMatrixTools::restrict_to_serial_sparse_matrix(
        laplace_matrix,
        sparsity_pattern,
        is_0,
        serial_sparse_matrix,
        serial_sparsity_pattern);
    else
      SparseMatrixTools::restrict_to_serial_sparse_matrix(
        laplace_matrix,
        sparsity_pattern,
        is_0,
        is_1,
        serial_sparse_matrix,
        serial_sparsity_pattern);

    FullMatrix<double> serial_sparse_matrix_full;
    serial_sparse_matrix_full.copy_from(serial_sparse_matrix);
    serial_sparse_matrix_full.print_formatted(deallog.get_file_stream(),
                                              2,
                                              false,
                                              8);
  };

  test_restrict(dof_handler.locally_owned_dofs(), {});
  test_restrict(DoFTools::extract_locally_active_dofs(dof_handler), {});
  test_restrict(dof_handler.locally_owned_dofs(),
                DoFTools::extract_locally_active_dofs(dof_handler));
}

#include "../tests.h"

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  // SparseMatrix -> SparseMatrix
  if (Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1)
    test<2,
         SparseMatrix<double>,
         SparsityPattern,
         SparseMatrix<double>,
         SparsityPattern>();

  // TrilinosWrappers::SparseMatrix -> SparseMatrix
  test<2,
       TrilinosWrappers::SparseMatrix,
       TrilinosWrappers::SparsityPattern,
       SparseMatrix<double>,
       SparsityPattern>();

  // TrilinosWrappers::SparseMatrix -> TrilinosWrappers::SparseMatrix
  test<2,
       TrilinosWrappers::SparseMatrix,
       TrilinosWrappers::SparsityPattern,
       TrilinosWrappers::SparseMatrix,
       TrilinosWrappers::SparsityPattern>();
}
