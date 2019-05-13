// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


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
