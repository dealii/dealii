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


// test setting some elements using a non-const chunk matrix iterator and
// operator+=, and reading them back through the matrix itself

#include <deal.II/lac/chunk_sparse_matrix.h>

#include "../tests.h"


void
test(const unsigned int chunk_size)
{
  deallog << "Chunk size: " << chunk_size << std::endl;
  ChunkSparsityPattern sp(5, 5, 3, chunk_size);
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j < 5; ++j)
      if (((i + 2 * j + 1) % 3 == 0) || (i == j))
        sp.add(i, j);
  sp.compress();

  ChunkSparseMatrix<double> m(sp);
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j < 5; ++j)
      if (((i + 2 * j + 1) % 3 == 0) || (i == j))
        m.set(i, j, 1.);

  ChunkSparseMatrix<double>::iterator i = m.begin();
  for (; i != m.end(); ++i)
    i->value() += i->row() * i->column();

  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j < 5; ++j)
      if (((i + 2 * j + 1) % 3 == 0) || (i == j))
        {
          deallog << i << ' ' << j << ' ' << m.el(i, j) << std::endl;
          Assert(std::fabs(m.el(i, j) - (1. + i * j)) < 1e-14,
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
      test(1);
      test(2);
      test(4);
      test(5);
      test(7);
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
