// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// FullMatrix::copy_from could not be compiled if we copied from a
// sparse matrix. make sure this now works

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>

#include "../tests.h"


int
main()
{
  initlog();

  const unsigned int N = 4;
  FullMatrix<double> f(N, N);

  SparsityPattern s(N, N, N);
  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < N; ++j)
      s.add(i, j);
  s.compress();

  SparseMatrix<double> sm(s);
  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < N; ++j)
      sm.set(i, j, i * j);

  f.copy_from(sm);

  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < N; ++j)
      {
        deallog << i << ' ' << j << ' ' << f(i, j) << std::endl;
        AssertThrow(f(i, j) == sm(i, j), ExcInternalError());
      }

  deallog << "OK" << std::endl;
}
