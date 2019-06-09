// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

// test read/write via RowAccessor
// similar to bcsr_09 but use a different sparsity to expose the bug
// that advance() did not "jump over" completely empty columns
// in sparsity below we have it as column 4

//     01   23  4   5  678
//      2    2  1   1   3
//  3   x           x        012
//  2        x          x    34
//  1   x           x   x    5
//  2   x    x               67

#include <deal.II/base/logstream.h>

#include <deal.II/lac/lapack_full_matrix.h>

#include <fstream>
#include <iostream>

#include "bcsr_helper.h"


using namespace dealii;

void
test()
{
  // number of blocks:
  const std::vector<unsigned int> row_blocks = {{3, 2, 1, 2}};
  const std::vector<unsigned int> col_blocks = {{2, 2, 1, 1, 3}};
  const unsigned int              M          = row_blocks.size();
  const unsigned int              N          = col_blocks.size();

  std::vector<dealii::types::global_dof_index> row_offset;
  std::vector<dealii::types::global_dof_index> col_offset;

  auto setup_offset = [](const std::vector<unsigned int> &             blocks,
                         std::vector<dealii::types::global_dof_index> &offset) {
    offset.resize(blocks.size() + 1, 0);
    std::partial_sum(blocks.begin(), blocks.end(), ++offset.begin());
  };

  setup_offset(row_blocks, row_offset);
  setup_offset(col_blocks, col_offset);

  deallog << "row blocks:";
  for (auto el : row_blocks)
    deallog << " " << el;
  deallog << std::endl;

  deallog << "col blocks:";
  for (auto el : col_blocks)
    deallog << " " << el;
  deallog << std::endl;

  deallog << "row offset:";
  for (auto el : row_offset)
    deallog << " " << el;
  deallog << std::endl;

  deallog << "col offset:";
  for (auto el : col_offset)
    deallog << " " << el;
  deallog << std::endl;

  DynamicSparsityPattern dsp(M, N);
  dsp.add(0, 0);
  dsp.add(0, 3);
  dsp.add(1, 1);
  dsp.add(1, 4);
  dsp.add(2, 0);
  dsp.add(2, 3);
  dsp.add(2, 4);
  dsp.add(3, 0);
  dsp.add(3, 1);

  std::shared_ptr<BlockIndices> rb = std::make_shared<BlockIndices>(row_blocks);
  std::shared_ptr<BlockIndices> cb = std::make_shared<BlockIndices>(col_blocks);

  auto bcsr_row_part =
    std::make_shared<dealii::Utilities::MPI::Partitioner>(rb->total_size());

  // setup matrices
  BlockCSRMatrix<double>        A;
  const BlockCSRMatrix<double> &A_const = A;
  A.reinit(dsp, rb, cb, bcsr_row_part);

  // setup
  init_bcsr(A);

  deallog << "m: " << A.m() << std::endl << "n: " << A.n() << std::endl;
  deallog << "initial:" << std::endl;
  const auto full_M =
    std::accumulate(row_blocks.begin(), row_blocks.end(), (unsigned int)0);
  const auto full_N =
    std::accumulate(col_blocks.begin(), col_blocks.end(), (unsigned int)0);

  deallog << "   ";
  for (unsigned int j = 0; j < full_N; ++j)
    deallog << "    " << j << "   ";
  deallog << std::endl;
  for (unsigned int i = 0; i < full_M; ++i)
    {
      deallog << i << "  ";
      for (unsigned int j = 0; j < full_N; ++j)
        deallog << " " << A_const.el(i, j);

      deallog << std::endl;
    }

  // now test:
  const std::vector<unsigned int> my_rows = {{1, 2, 3, 7}};

  bool iterate = true;
  auto read    = [&]() -> void {
    BlockCSRMatrixIterators::RowsAccessor<double, true> const_row_accessor(
      &A, my_rows);
    deallog << std::endl << "Reading: " << std::endl;
    while (iterate)
      {
        deallog << std::endl
                << "Column: " << const_row_accessor.current_column()
                << std::endl;
        for (auto ind : my_rows)
          deallog << ind << "  " << const_row_accessor.local_element(ind)
                  << std::endl;

        iterate = const_row_accessor.advance();
      }
  };

  read();

  iterate = true;
  {
    BlockCSRMatrixIterators::RowsAccessor<double, false> row_accessor(&A,
                                                                      my_rows);
    deallog << std::endl << "Setting to zero: " << std::endl;
    while (iterate)
      {
        deallog << "Column: " << row_accessor.current_column() << std::endl;
        for (auto ind : my_rows)
          row_accessor.set_local_element(ind, 0.);

        iterate = row_accessor.advance();
      }
  }

  iterate = true;
  read();

  deallog << "Ok" << std::endl;
}

int
main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  std::ofstream                            logfile("output");
  dealii::deallog.attach(logfile, /*do not print job id*/ false);
  dealii::deallog.depth_console(0);

  test();
}
