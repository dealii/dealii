// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this tests a failure in the design of the block sparse matrix iterators:
// falling off the end of the matrix does not yield the iterator provided by the
// end() function

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  // create a 1x2 block matrix with
  // non-quadratic blocks so as to make them
  // not specially store the diagonal
  BlockSparsityPattern bsp(1, 2);
  for (unsigned int j = 0; j < 2; ++j)
    bsp.block(0, j).reinit(3, 2, 1);
  bsp.collect_sizes();

  // leave row 0 of block 0,0 empty, but have
  // something in this row for block 0,1
  bsp.block(0, 0).add(1, 0);
  bsp.block(0, 1).add(0, 0);
  bsp.compress();

  BlockSparseMatrix<double> m(bsp);

  // get the start iterator. it should point
  // to the first element of the first row,
  // which happens to be in block 0,1
  BlockSparseMatrix<double>::const_iterator it = m.begin();

  deallog << it->row() << ' ' << it->column() << ' ' << it->block_row() << ' '
          << it->block_column() << std::endl;

  Assert(it->row() == 0, ExcInternalError());
  Assert(it->column() == 2, ExcInternalError());
  Assert(it->block_row() == 0, ExcInternalError());
  Assert(it->block_column() == 1, ExcInternalError());

  // now advance by two (the only two
  // elements of the matrix) and make sure
  // that we equal the end iterator
  ++it;
  ++it;
  Assert(it == m.end(), ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  try
    {
      test();
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
