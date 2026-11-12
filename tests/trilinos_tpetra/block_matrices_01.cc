// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2022 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Verify n_nonzero_elements and reinit() for Tpetra BlockSparseMatrix.

#include <deal.II/base/mpi.h>

#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/trilinos_tpetra_block_sparse_matrix.h>

#include "../tests.h"



template <typename Number, typename MemorySpace>
void
test()
{
  // create block sparsity pattern
  BlockDynamicSparsityPattern bdsp(2, 2);
  bdsp.block(0, 0).reinit(2, 2);
  bdsp.block(0, 1).reinit(2, 3);
  bdsp.block(1, 0).reinit(3, 2);
  bdsp.block(1, 1).reinit(3, 3);
  bdsp.collect_sizes();

  // add triangular sparsity pattern to each block
  // block (0, 0)
  for (unsigned int row = 0; row < 2; ++row)
    for (unsigned int col = 0; col <= row; ++col)
      bdsp.add(row, col);
  // block (1, 1)
  for (unsigned int row = 2; row < 5; ++row)
    for (unsigned int col = 2; col <= row; ++col)
      bdsp.add(row, col);

  bdsp.compress();
  AssertThrow(bdsp.n_nonzero_elements() == 9, ExcInternalError());

  // create block sparse matrix
  LinearAlgebra::TpetraWrappers::BlockSparseMatrix<Number, MemorySpace> tbsm;
  tbsm.reinit(bdsp);
  AssertThrow(tbsm.n_nonzero_elements() == 9, ExcInternalError());



  // create partitioning
  std::vector<dealii::IndexSet> partitioning(2);
  partitioning[0] = dealii::IndexSet(2); // 2 dofs
  partitioning[1] = dealii::IndexSet(3); // 3 dofs

  partitioning[0].add_range(0, 2);
  partitioning[1].add_range(0, 3);

  for (auto &set : partitioning)
    set.compress();

  // create block sparse matrix using another reinit() function
  LinearAlgebra::TpetraWrappers::BlockSparseMatrix<Number, MemorySpace> tbsm2;
  tbsm2.reinit(partitioning, bdsp, MPI_COMM_SELF, false);

  AssertThrow(tbsm.n_nonzero_elements() == 9, ExcInternalError());
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

#ifdef DEAL_II_TRILINOS_WITH_TPETRA_INST_DOUBLE
  test<double, MemorySpace::Host>();
  test<double, MemorySpace::Default>();
#endif

#ifdef DEAL_II_TRILINOS_WITH_TPETRA_INST_FLOAT
  test<float, MemorySpace::Host>();
  test<float, MemorySpace::Default>();
#endif

  deallog << "OK" << std::endl;

  return 0;
}
