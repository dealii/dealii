// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Check TpetraWrappers::SparseMatrix::add() with elide_zero_values == true.

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

#include "../tests.h"


template <typename MemorySpace>
void
test()
{
  const unsigned int n = 3;

  DynamicSparsityPattern dsp(n, n);
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
      dsp.add(i, j);

  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace> matrix;
  matrix.reinit(dsp);

  const types::global_dof_index columns[n] = {0, 1, 2};
  const double values[n][n]                = {{1, 2, 3}, {4, 0, 5}, {0, 0, 0}};
  for (unsigned int i = 0; i < n; ++i)
    matrix.add(i, n, columns, values[i], /* elide_zero_values = */ true);
  matrix.compress(VectorOperation::add);

  for (unsigned int i = 0; i < n; ++i)
    deallog << matrix(i, 0) << ' ' << matrix(i, 1) << ' ' << matrix(i, 2)
            << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  initlog();

  test<MemorySpace::Host>();
  test<MemorySpace::Default>();
}
