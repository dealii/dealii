// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// test LinearAlgebra::TpetraWrappers::SparseMatrix<double> manual assembly on
// device

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_sparsity_pattern.h>

#include <iostream>

#include "../tests.h"


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  initlog();

  int N = 50;

  SparsityPattern sparsity(N, N, 5);
  for (int i = 0; i < N; ++i)
    {
      if (i > 0)
        sparsity.add(i, i - 1);
      sparsity.add(i, i);
      if (i < N - 1)
        sparsity.add(i, i + 1);
    }
  sparsity.compress();
  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default>
    matrix;
  matrix.reinit(sparsity);

  {
    auto localMatrix = matrix.trilinos_matrix().getLocalMatrixDevice();

    matrix.trilinos_matrix().resumeFill();


    Kokkos::parallel_for(
      Kokkos::RangePolicy<MemorySpace::Default::kokkos_space::execution_space>(
        1, N - 1),
      KOKKOS_LAMBDA(const size_t &element_gidx) {
        int local_row_id = element_gidx;

        int col_idx[3] = {local_row_id, local_row_id - 1, local_row_id + 1};

        int    n_values = 3;
        double values[] = {0.11, 0.11, 0.11};

        localMatrix.sumIntoValues(
          local_row_id, col_idx, n_values, values, true, true);
        n_values = 1;
        localMatrix.sumIntoValues(
          local_row_id - 1, col_idx, n_values, values, true, true);
        localMatrix.sumIntoValues(
          local_row_id + 1, col_idx, n_values, values, true, true);
      });

    matrix.trilinos_matrix().fillComplete();
  }

  deallog << "matrix:" << std::endl;
  matrix.print(deallog.get_file_stream());
}
