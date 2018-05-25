// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// this tests a failure in the design of the block sparse matrix iterators:
// falling off the end of the matrix does not yield the iterator provided by the
// end() function

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  BlockSparsityPattern bsp(2, 2);
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      bsp.block(i, j).reinit(1, 1, 1);
  bsp.collect_sizes();
  bsp.compress();

  BlockSparseMatrix<double> m(bsp);

  // advance it to the end of the matrix
  BlockSparseMatrix<double>::const_iterator it = m.begin();
  for (unsigned int i = 0; i < 4; ++i)
    ++it;

  // now also get an end iterator
  BlockSparseMatrix<double>::const_iterator it2 = m.end();

  // make sure that the two of them match
  AssertThrow(it == it2, ExcInternalError());

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
  catch (std::exception &exc)
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
