// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Tests LAPACKFullMatrix::operator*= and operator/=

#include <deal.II/lac/lapack_full_matrix.h>

#include <iostream>

#include "../tests.h"



void
test()
{
  const unsigned int       m = 7;
  const unsigned int       n = 11;
  LAPACKFullMatrix<double> A(m, n);
  for (unsigned int i = 0; i < m; ++i)
    for (unsigned int j = 0; j < n; ++j)
      A(i, j) = random_value<double>();

  LAPACKFullMatrix<double> A_check(A);

  const double factor = 2.345;

  // multiply by factor
  A *= factor;

  // test multiplication
  for (unsigned int i = 0; i < m; ++i)
    for (unsigned int j = 0; j < n; ++j)
      AssertThrow(std::abs(A(i, j) - factor * A_check(i, j)) < 1.e-12,
                  ExcInternalError());

  // divide by factor
  A /= factor;

  // test division
  for (unsigned int i = 0; i < m; ++i)
    for (unsigned int j = 0; j < n; ++j)
      AssertThrow(std::abs(A(i, j) - A_check(i, j)) < 1.e-12,
                  ExcInternalError());

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  test();
}
