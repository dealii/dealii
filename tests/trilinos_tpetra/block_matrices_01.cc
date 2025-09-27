// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Verify n_nonzero_elements and reinit() for Tpetra BlockSparseMatrix.

#include <deal.II/base/mpi.h>

#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>

#include "../tests.h"



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
  deallog << "nonzeros BlockSparsityPattern: " << bdsp.n_nonzero_elements()
          << std::endl;

  // create block sparse matrix
  TrilinosWrappers::BlockSparseMatrix tbsm;
  tbsm.reinit(bdsp);
  deallog << "nonzeros BlockSparseMatrix: " << tbsm.n_nonzero_elements()
          << std::endl;



  // create partitioning
  std::vector<dealii::IndexSet> partitioning(2);
  partitioning[0] = dealii::IndexSet(2); // 2 dofs
  partitioning[1] = dealii::IndexSet(3); // 3 dofs

  partitioning[0].add_range(0, 2);
  partitioning[1].add_range(0, 3);

  for (auto &set : partitioning)
    set.compress();

  // create block sparse matrix using another reinit() function
  TrilinosWrappers::BlockSparseMatrix tbsm2;
  tbsm2.reinit(partitioning, bdsp, MPI_COMM_SELF, false);

  deallog << "nonzeros BlockSparseMatrix: " << tbsm.n_nonzero_elements()
          << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  test();

  return 0;
}
