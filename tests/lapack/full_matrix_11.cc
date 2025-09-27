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


// Tests LAPACKFullMatrix::TmTmult

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>

#include "../tests.h"



void
test()
{
  const unsigned int       m = 2;
  const unsigned int       n = 3;
  const unsigned int       k = 4;
  FullMatrix<double>       A(k, m), B(n, k), C(m, n), OC(m, n);
  LAPACKFullMatrix<double> AL(k, m), BL(n, k), CL(m, n);
  for (unsigned int i = 0; i < k; ++i)
    for (unsigned int j = 0; j < m; ++j)
      A(i, j) = AL(i, j) = random_value<double>();
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < k; ++j)
      B(i, j) = BL(i, j) = random_value<double>();

  A.TmTmult(C, B);
  AL.TmTmult(CL, BL);
  AL.TmTmult(OC, BL);
  for (unsigned int i = 0; i < m; ++i)
    for (unsigned int j = 0; j < n; ++j)
      {
        Assert(std::abs(C(i, j) - CL(i, j)) < 1e-13, ExcInternalError());
        Assert(std::abs(C(i, j) - OC(i, j)) < 1e-13, ExcInternalError());
      }

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  test();
}
