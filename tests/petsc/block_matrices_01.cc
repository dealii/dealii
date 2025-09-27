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
  Assert(tmp2.petsc_matrix() == pbsm.petsc_matrix(), ExcInternalError());
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

  // Extract the PETSc MATNEST, create an array of PETSc matrices and assign
  // to a new BlockSparseMatrix
  std::array<Mat, 2> arrayRows0 = {
    {tmp2.block(0, 0).petsc_matrix(), tmp2.block(0, 1).petsc_matrix()}};
  std::array<Mat, 2> arrayRows1 = {
    {tmp2.block(1, 0).petsc_matrix(), tmp2.block(1, 1).petsc_matrix()}};
  std::array<std::array<Mat, 2>, 2> arrayMat = {{arrayRows0, arrayRows1}};

  PETScWrappers::MPI::BlockSparseMatrix tmp3(arrayMat);
  Assert(tmp3.n_block_rows() == pbsm.n_block_rows(), ExcInternalError());
  Assert(tmp3.n_block_cols() == pbsm.n_block_cols(), ExcInternalError());
  Assert(tmp3.m() == pbsm.m(), ExcInternalError());
  Assert(tmp3.n() == pbsm.n(), ExcInternalError());
  for (unsigned int blr = 0; blr < 2; ++blr)
    {
      for (unsigned int blc = 0; blc < 2; ++blc)
        {
          Assert(tmp3.block(blr, blc).m() == pbsm.block(blr, blc).m(),
                 ExcInternalError());
          Assert(tmp3.block(blr, blc).n() == pbsm.block(blr, blc).n(),
                 ExcInternalError());
          Assert(tmp3.block(blr, blc).petsc_matrix() ==
                   pbsm.block(blr, blc).petsc_matrix(),
                 ExcInternalError());
        }
    }

  // Now pass empty blocks
  std::array<Mat, 2> arrayRows0Empty = {
    {nullptr, tmp2.block(0, 1).petsc_matrix()}};
  std::array<Mat, 2> arrayRows1Empty = {
    {tmp2.block(1, 0).petsc_matrix(), nullptr}};
  std::array<std::array<Mat, 2>, 2> arrayMatEmpty = {
    {arrayRows0Empty, arrayRows1Empty}};

  PETScWrappers::MPI::BlockSparseMatrix tmp4(arrayMatEmpty);
  Assert(tmp4.n_block_rows() == pbsm.n_block_rows(), ExcInternalError());
  Assert(tmp4.n_block_cols() == pbsm.n_block_cols(), ExcInternalError());
  Assert(tmp4.m() == pbsm.m(), ExcInternalError());
  Assert(tmp4.n() == pbsm.n(), ExcInternalError());
  Assert(tmp4.block(0, 1).m() == pbsm.block(0, 1).m(), ExcInternalError());
  Assert(tmp4.block(0, 1).n() == pbsm.block(0, 1).n(), ExcInternalError());
  Assert(tmp4.block(0, 1).petsc_matrix() == pbsm.block(0, 1).petsc_matrix(),
         ExcInternalError());
  Assert(tmp4.block(1, 0).m() == pbsm.block(1, 0).m(), ExcInternalError());
  Assert(tmp4.block(1, 0).n() == pbsm.block(1, 0).n(), ExcInternalError());
  Assert(tmp4.block(1, 0).petsc_matrix() == pbsm.block(1, 0).petsc_matrix(),
         ExcInternalError());

  // Check the rectangular cases
  std::array<std::array<Mat, 2>, 1> arrayMatRow = {{arrayRows0}};
  std::array<std::array<Mat, 1>, 2> arrayMatCol = {
    {{{arrayMat[0][0]}}, {{arrayMat[1][0]}}}};

  PETScWrappers::MPI::BlockSparseMatrix tmp5(arrayMatRow);
  Assert(tmp5.n_block_cols() == pbsm.n_block_cols(), ExcInternalError());
  Assert(tmp5.n() == pbsm.n(), ExcInternalError());
  for (unsigned int blc = 0; blc < 2; ++blc)
    {
      Assert(tmp5.block(0, blc).m() == pbsm.block(0, blc).m(),
             ExcInternalError());
      Assert(tmp5.block(0, blc).n() == pbsm.block(0, blc).n(),
             ExcInternalError());
      Assert(tmp5.block(0, blc).petsc_matrix() ==
               pbsm.block(0, blc).petsc_matrix(),
             ExcInternalError());
    }

  PETScWrappers::MPI::BlockSparseMatrix tmp6(arrayMatCol);
  Assert(tmp6.n_block_rows() == pbsm.n_block_rows(), ExcInternalError());
  Assert(tmp6.m() == pbsm.m(), ExcInternalError());
  for (unsigned int blr = 0; blr < 2; ++blr)
    {
      Assert(tmp6.block(blr, 0).m() == pbsm.block(blr, 0).m(),
             ExcInternalError());
      Assert(tmp6.block(blr, 0).n() == pbsm.block(blr, 0).n(),
             ExcInternalError());
      Assert(tmp6.block(blr, 0).petsc_matrix() ==
               pbsm.block(blr, 0).petsc_matrix(),
             ExcInternalError());
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
