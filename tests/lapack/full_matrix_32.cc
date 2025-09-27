// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Tests LAPACKFullMatrix::scale_rows()

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>
#include <tuple>

#include "../tests.h"

#include "create_matrix.h"

DeclException5(ExcEl,
               int,
               int,
               double,
               double,
               double,
               << "Error in element (" << arg1 << ',' << arg2 << "): " << arg3
               << " != " << arg4 << " delta=" << arg5);

template <typename NumberType>
void
test(const unsigned int n, const unsigned int k, const NumberType eps)
{
  deallog << n << ' ' << k << ' ' << std::endl;
  FullMatrix<NumberType>       A(k, n), C(k, n), D(k, k);
  LAPACKFullMatrix<NumberType> AL(k, n);
  Vector<NumberType>           DL(k);

  create_random(AL);
  create_random(DL);

  A = AL;
  D = NumberType();
  for (unsigned int i = 0; i < k; ++i)
    D(i, i) = DL(i);

  // C = D*A
  D.mmult(C, A);

  // A = D * A
  AL.scale_rows(DL);

  for (unsigned int i = 0; i < k; ++i)
    for (unsigned int j = 0; j < n; ++j)
      AssertThrow(std::abs(C(i, j) - AL(i, j)) < eps * std::abs(AL(i, j)),
                  ExcEl(i, j, C(i, j), AL(i, j), C(i, j) - AL(i, j)));

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  const std::vector<std::array<unsigned int, 2>> sizes = {{{3, 3}},
                                                          {{7, 7}},
                                                          {{51, 51}},
                                                          {{320, 320}},
                                                          {{3, 9}},
                                                          {{9, 7}},
                                                          {{5, 17}},
                                                          {{320, 121}}};

  deallog.push("double");
  for (auto el : sizes)
    test<double>(el[0], el[1], 1e-13);
  deallog.pop();

  deallog.push("float");
  for (auto el : sizes)
    test<float>(el[0], el[1], 1e-5);
  deallog.pop();
}
