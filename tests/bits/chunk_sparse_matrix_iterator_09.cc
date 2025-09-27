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



// this test is sparse_matrix_iterator_09 for a ChunkSparseMatrix

#include <deal.II/lac/chunk_sparse_matrix.h>

#include "../tests.h"


void
test(const unsigned int chunk_size)
{
  deallog << "Chunk size: " << chunk_size << std::endl;

  // create a sparsity pattern with totally
  // empty lines (not even diagonals, since
  // not quadratic)
  ChunkSparsityPattern sparsity(4, 5, 1, chunk_size);
  sparsity.add(1, 1);
  sparsity.add(3, 1);
  sparsity.compress();

  // attach a sparse matrix to it
  ChunkSparseMatrix<double> A(sparsity);

  // and loop over the elements of it
  for (ChunkSparseMatrix<double>::const_iterator k = A.begin(); k != A.end();
       ++k)
    deallog << k->row() << ' ' << k->column() << ' ' << k->value() << std::endl;
}



int
main()
{
  initlog();

  try
    {
      test(1);
      test(3);
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
