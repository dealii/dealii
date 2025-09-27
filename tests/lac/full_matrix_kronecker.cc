// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test that an assertion is thrown if output argument of FullMatrix::mmult() is
// also an input argument

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog << std::fixed;

  dealii::FullMatrix<double> A(2, 2);
  dealii::FullMatrix<double> B(2, 2);
  dealii::FullMatrix<double> C(4, 4);

  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  B(0, 0) = 5;
  B(0, 1) = 6;
  B(1, 0) = 7;
  B(1, 1) = 8;

  C.kronecker_product(A, B);


  std::stringstream sstring;
  sstring << "Result"
          << "\n";
  C.print_formatted(sstring, 4, true, 10, "0.");
  deallog << sstring.str() << std::endl;

  return 0;
}
