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


// Tests LAPACKFullMatrix::set(i,j,a)

#include <deal.II/lac/lapack_full_matrix.h>

#include <iostream>
#include <tuple>

#include "../tests.h"

template <typename NumberType>
void
test(const unsigned int n = 3, const unsigned int k = 6)
{
  LAPACKFullMatrix<NumberType> A(n, k);
  A.set(0, 1, 2.);
  AssertThrow(A(0, 1) == 2., ExcInternalError());

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  test<double>();
}
