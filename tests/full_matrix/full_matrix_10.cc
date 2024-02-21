// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check method mTmult of FullMatrix on larger size than full_matrix_02 where
// we interface to the external BLAS

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"


template <typename Number>
void
test()
{
  FullMatrix<Number> A(2, 76), B(3, 76), C(2, 3), D(2, 3);
  for (unsigned int i = 0; i < A.m(); ++i)
    for (unsigned int j = 0; j < A.n(); ++j)
      A(i, j) = random_value<double>();
  for (unsigned int i = 0; i < B.m(); ++i)
    for (unsigned int j = 0; j < B.n(); ++j)
      B(i, j) = random_value<double>();

  A.mTmult(C, B); // C = A * B^T

  for (unsigned int i = 0; i < A.m(); ++i)
    for (unsigned int j = 0; j < B.m(); ++j)
      for (unsigned int k = 0; k < A.n(); ++k)
        D(i, j) += A(i, k) * B(j, k);
  C.add(-1., D);

  deallog << "Difference: "
          << filter_out_small_numbers(C.l1_norm(),
                                      std::numeric_limits<Number>::epsilon() *
                                        100.)
          << std::endl;
}

int
main()
{
  initlog();

  test<double>();
  test<float>();
  test<long double>();
}
