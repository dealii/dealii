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



// check ChunkSparseMatrix::vmult

#include <deal.II/lac/chunk_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <vector>

#include "../tests.h"


void
test(const unsigned int chunk_size, Vector<double> &v, Vector<double> &w)
{
  // set some entries in the
  // matrix. actually, set them all
  ChunkSparsityPattern sp(v.size(), v.size(), v.size(), chunk_size);
  for (unsigned int i = 0; i < v.size(); ++i)
    for (unsigned int j = 0; j < v.size(); ++j)
      sp.add(i, j);
  sp.compress();

  // then create a matrix from that
  ChunkSparseMatrix<double> m(sp);
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.n(); ++j)
      m.set(i, j, i + 2 * j);

  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = i;

  v.compress();
  w.compress();

  // w:=Mv
  m.vmult(w, v);

  // make sure we get the expected result
  for (unsigned int i = 0; i < v.size(); ++i)
    {
      Assert(v(i) == i, ExcInternalError());

      double result = 0;
      for (unsigned int j = 0; j < m.n(); ++j)
        result += (i + 2 * j) * j;
      Assert(w(i) == result, ExcInternalError());
    }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  try
    {
      const unsigned int chunk_sizes[] = {1, 2, 4, 7, 11};
      for (unsigned int i = 0; i < sizeof(chunk_sizes) / sizeof(chunk_sizes[0]);
           ++i)
        {
          Vector<double> v(100);
          Vector<double> w(100);
          test(chunk_sizes[i], v, w);
        }
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
