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



// set a few elements in a chunk sparse matrix and test for iterator
// inequality

#include <deal.II/lac/chunk_sparse_matrix.h>

#include "../tests.h"


void
test(const unsigned int chunk_size)
{
  deallog << "Chunk size = " << chunk_size << std::endl;

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

  // then extract the elements (note that
  // some may be zero or even outside the
  // matrix
  AssertDimension(m.end() - m.begin(), m.n_nonzero_elements());
  for (unsigned int i = 0; i < m.m(); ++i)
    {
      deallog << "row " << i << ": ";
      AssertDimension(m.end(i) - m.begin(i),
                      m.get_sparsity_pattern().row_length(i));
      for (ChunkSparseMatrix<double>::const_iterator it = m.begin(i);
           it != m.end(i);
           ++it)
        {
          deallog << " done " << (it - m.begin(i)) << ", left "
                  << (it - m.end(i));
        }
      deallog << std::endl;
    }
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
