// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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


// Verify n_nonzero_elements for PETSc BlockSparseMatrix.

#include <deal.II/base/mpi.h>

#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/petsc_block_sparse_matrix.h>

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
  std::vector<IndexSet> rows(2);
  rows[0].set_size(2);
  rows[0].add_range(0, 2);
  rows[1].set_size(3);
  rows[1].add_range(0, 3);

  std::vector<IndexSet> cols(2);
  cols[0].set_size(2);
  cols[0].add_range(0, 2);
  cols[1].set_size(3);
  cols[1].add_range(0, 3);

  PETScWrappers::MPI::BlockSparseMatrix pbsm;
  pbsm.reinit(rows, cols, bdsp, MPI_COMM_WORLD);
  deallog << "nonzeros BlockSparseMatrix: " << pbsm.n_nonzero_elements()
          << std::endl;

  // Extract the PETSc MATNEST and use print from PETScWrappers::MatrixBase
  PETScWrappers::MatrixBase tmp(pbsm.petsc_matrix());
  tmp.print(deallog.get_file_stream());

  // Extract the PETSc MATNEST and assign to a new BlockSparseMatrix
  PETScWrappers::MPI::BlockSparseMatrix tmp2(pbsm.petsc_matrix());
  Assert(tmp2.n_block_rows() == pbsm.n_block_rows(), ExcInternalError());
  Assert(tmp2.n_block_cols() == pbsm.n_block_cols(), ExcInternalError());
  Assert(tmp2.m() == pbsm.m(), ExcInternalError());
  Assert(tmp2.n() == pbsm.n(), ExcInternalError());
  for (unsigned int blr = 0; blr < 2; ++blr)
    {
      for (unsigned int blc = 0; blc < 2; ++blc)
        {
          Assert(tmp2.block(blr, blc).m() == pbsm.block(blr, blc).m(),
                 ExcInternalError());
          Assert(tmp2.block(blr, blc).n() == pbsm.block(blr, blc).n(),
                 ExcInternalError());
          Assert(tmp2.block(blr, blc).petsc_matrix() ==
                   pbsm.block(blr, blc).petsc_matrix(),
                 ExcInternalError());
        }
    }
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  test();

  return 0;
}
