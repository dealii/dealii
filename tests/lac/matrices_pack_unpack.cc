// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Copying from a matrix into std::vector and vice versa using iterators

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"

int
main()
{
  initlog();

  dealii::FullMatrix<double> A(3, 4);
  std::iota(A.begin(), A.end(), 1.0);

  deallog << "A =" << std::endl;
  A.print(deallog);

  std::vector<double> v(A.n_elements());
  std::copy(A.begin(), A.end(), v.begin());
  deallog << "v = " << std::endl;
  deallog << v << std::endl;

  dealii::FullMatrix<double> B(A.m(), A.n());
  std::copy(v.cbegin(), v.cend(), B.begin());
  deallog << "B =" << std::endl;
  A.print(deallog);

  return 0;
}
