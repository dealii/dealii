// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Tests reinitialization of square and rectangle LAPACKFullMatrix

#include <deal.II/lac/lapack_full_matrix.h>

#include <iostream>

#include "../tests.h"


void
test(const unsigned int size, const bool reinit_square)
{
  // this test can not currently work with matrices smaller than
  // 1\times2.
  AssertThrow(size > 2, ExcInternalError());

  // initialise a first matrix with the standard constructor and fill
  // it with some numbers
  LAPACKFullMatrix<double> M(size, size);

  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      M(i, j) = i + 2. * j;

  // initialise a second matrix with the standard constructor and fill
  // it with some numbers
  LAPACKFullMatrix<double> N(size + 2, size - 2);

  for (unsigned int i = 0; i < N.m(); ++i)
    for (unsigned int j = 0; j < N.n(); ++j)
      N(i, j) = i + 2. * j;

  // clearly, this should be the case
  AssertThrow(N.m() != M.m(), ExcInternalError());
  AssertThrow(N.n() != M.n(), ExcInternalError());

  // if reinit_square is true, reinitialise the rectangle matrix to a
  // square matrix (use reinit (const unsigned int))
  if (reinit_square)
    {
      // reinitialise the matrix and fill it with some numbers
      N.reinit(size);

      for (unsigned int i = 0; i < N.m(); ++i)
        for (unsigned int j = 0; j < N.n(); ++j)
          N(i, j) = i + 2. * j;
    }

  // otherwise reinitialise the rectangle matrix to a square one (use
  // reinit (const unsigned int, const unsigned int))
  else
    {
      // reinitialise the matrix and fill it with some numbers
      M.reinit(size + 2, size - 2);

      for (unsigned int i = 0; i < M.m(); ++i)
        for (unsigned int j = 0; j < M.n(); ++j)
          M(i, j) = i + 2. * j;
    }

  // and now this should be true
  AssertThrow(N.m() == M.m(), ExcInternalError());
  AssertThrow(N.n() == M.n(), ExcInternalError());

  // in fact, this should be true too, so check
  for (unsigned int i = 0; i < M.m(); ++i)
    for (unsigned int j = 0; j < M.n(); ++j)
      AssertThrow(M(i, j) == N(i, j), ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  // Test square matrix initialization
  test(4, true);
  test(5, true);
  test(6, true);

  // Test rectangle matrix initialization
  test(4, false);
  test(5, false);
  test(6, false);
}
