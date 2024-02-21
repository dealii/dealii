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



// check setting elements in a sparse matrix using ChunkSparseMatrix::set().
// make sure they are correct, and make sure that for the nonexisting entries
// ChunkSparseMatrix::el() returns zero and operator() throws an exception

#include <deal.II/lac/chunk_sparse_matrix.h>

#include "../tests.h"


void
test(const unsigned int chunk_size)
{
  ChunkSparsityPattern sp(5, 5, 3, chunk_size);
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j < 5; ++j)
      if ((i + 2 * j + 1) % 3 == 0)
        sp.add(i, j);
  sp.compress();

  ChunkSparseMatrix<double> m(sp);

  // first set a few entries
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.n(); ++j)
      if ((i + 2 * j + 1) % 3 == 0)
        m.set(i, j, i * j * .5 + .5);

  // then make sure we retrieve the same ones
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.n(); ++j)
      if ((i + 2 * j + 1) % 3 == 0)
        {
          Assert(m(i, j) == i * j * .5 + .5, ExcInternalError());
          Assert(m.el(i, j) == i * j * .5 + .5, ExcInternalError());
        }
      else
        {
          // reading elements not in the
          // sparsity pattern should return
          // zero
          const double x = m.el(i, j);
          Assert(x == 0, ExcInternalError());

          // if this is a sparsity_pattern
          // with chunk_size==1, then we need
          // to get an exception if we access
          // any other element. if
          // chunk_size>1, then this isn't
          // necessarily true. because the
          // matrix is square, we may also
          // always access the diagonal
          // element
          bool exc_thrown = false;
          try
            {
              m(i, j);
            }
          catch (const std::exception &)
            {
              exc_thrown = true;
            }
          Assert((exc_thrown == true) || (chunk_size > 1) || (i == j),
                 ExcInternalError());
        }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  try
    {
      const unsigned int chunk_sizes[] = {1, 2, 4, 5, 7};
      for (unsigned int i = 0; i < sizeof(chunk_sizes) / sizeof(chunk_sizes[0]);
           ++i)
        test(chunk_sizes[i]);
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
