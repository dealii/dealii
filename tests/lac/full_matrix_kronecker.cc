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
  dealii::FullMatrix<double> C(2, 2);
  dealii::FullMatrix<double> result1(4, 4);
  dealii::FullMatrix<double> result2(4, 4);

  dealii::FullMatrix<double> result3;
  dealii::FullMatrix<double> result3_ref;

  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;

  B(0, 0) = 5;
  B(0, 1) = 6;
  B(1, 0) = 7;
  B(1, 1) = 8;

  C(0, 0) = 9;
  C(0, 1) = 10;
  C(1, 0) = 11;
  C(1, 1) = 12;


  result1.kronecker_product(A, B);

  for (unsigned int i = 0; i < result2.m(); ++i)
    for (unsigned int j = 0; j < result2.n(); ++j)
      result2(i, j) = 1.0;

  result2.kronecker_product(A, B, true);

  for (unsigned int i = 0; i < result2.m(); ++i)
    for (unsigned int j = 0; j < result2.n(); ++j)
      {
        result2(i, j) -= 1.0;
        {
          if (std::abs(result1(i, j) - result2(i, j)) > 1e-12)
            {
              deallog << "Kronecker product with adding=true failed at (" << i
                      << "," << j << "): " << result1(i, j) << " vs "
                      << result2(i, j) << std::endl;
            }
        }
      }


  std::stringstream sstring;
  sstring << "Result"
          << "\n";
  result2.print_formatted(sstring, 4, true, 10, "0.");
  deallog << sstring.str() << std::endl;

  result3.kronecker_product(A, B, C);
  result3_ref.kronecker_product(result1, C);

  // Check whether result3 and result3_ref are identical
  if (result3.m() != result3_ref.m() || result3.n() != result3_ref.n())
    {
      deallog << "result3 and result3_ref have different sizes: ("
              << result3.m() << "," << result3.n() << ") vs ("
              << result3_ref.m() << "," << result3_ref.n() << ")" << std::endl;
    }
  else
    {
      for (unsigned int i = 0; i < result3.m(); ++i)
        for (unsigned int j = 0; j < result3.n(); ++j)
          {
            const double diff = std::abs(result3(i, j) - result3_ref(i, j));
            if (diff > 1e-12)
              {
                deallog << "result3 differs from result3_ref at (" << i << ","
                        << j << "): " << result3(i, j) << " vs "
                        << result3_ref(i, j) << " (diff=" << diff << ")"
                        << std::endl;
              }
          }
    }


  return 0;
}
