// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include "../tests.h"

// A.fill() produced an ExcIndexRange(r,0,m()) exception with
// the additional Information: Index 6 is not in [0,3[.
// Bug reported by Florian Prill

int
main()
{
  initlog();

  // matrix sizes
  const unsigned int m = 3;
  const unsigned int n = 10;

  LAPACKFullMatrix<double> A(n);
  FullMatrix<double>       C(m);
  // fill some entries:
  C(0, 0)         = 1.0;
  C(m - 1, m - 1) = 1.0;
  // insert C into A's middle:
  A.fill(C, 3, 3, 0, 0);
  // check some values
  Assert(A(3, 3) == 1, ExcInternalError());
  Assert(A(5, 5) == 1, ExcInternalError());

  deallog << "OK" << std::endl;
}
