// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check method Tmmult of FullMatrix, symmetric case

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"

const double entries_A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
const double compare[9]   = {66, 78, 90, 78, 93, 108, 90, 108, 126};

int
main()
{
  initlog();
  deallog << std::fixed;
  deallog << std::setprecision(3);

  FullMatrix<double> A(3, 3, entries_A);
  FullMatrix<double> C(3, 3);
  FullMatrix<double> D(3, 3, compare);

  // compute C= A^T*A
  A.Tmmult(C, A);

  C.add(-1., D);
  Assert(C.frobenius_norm() < 1e-12, ExcInternalError());

  deallog << "OK" << std::endl;
}
