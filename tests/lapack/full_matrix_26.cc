// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Tests grow_or_shrink() of square and rectangle LAPACKFullMatrix

#include <deal.II/lac/lapack_full_matrix.h>

#include <iostream>

#include "../tests.h"


void
test(const unsigned int size)
{
  AssertThrow(size > 2, ExcInternalError());
  const unsigned int smaller = size - 2;
  const unsigned int larger  = size + 3;

  // initialise a first matrix with the standard constructor and fill
  // it with some numbers
  LAPACKFullMatrix<double> M(size, size);

  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      M(i, j) = i + 2. * j;

  LAPACKFullMatrix<double> M1(M), M2(M);

  M1.grow_or_shrink(smaller);
  for (unsigned int i = 0; i < smaller; ++i)
    for (unsigned int j = 0; j < smaller; ++j)
      AssertThrow(M1(i, j) == M(i, j), ExcInternalError());

  M2.grow_or_shrink(larger);
  for (unsigned int i = 0; i < larger; ++i)
    for (unsigned int j = 0; j < larger; ++j)
      {
        if (i < size && j < size)
          {
            AssertThrow(M2(i, j) == M(i, j), ExcInternalError());
          }
        else
          {
            AssertThrow(M2(i, j) == 0., ExcInternalError());
          }
      }

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  test(4);
  test(7);
  test(11);
  test(31);
}
