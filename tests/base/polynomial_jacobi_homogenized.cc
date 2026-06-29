// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// check jacobi_polynomial_homogenized_value by comparing it to its definition

#include <deal.II/base/polynomial.h>

#include "../tests.h"

using namespace Polynomials;


void
test()
{
  for (int alpha = 0; alpha < 3; ++alpha)
    for (int beta = 0; beta < 3; ++beta)
      for (unsigned int degree = 0; degree < 6; ++degree)
        for (const double x : {0.3, 1.0 / 3.0, 0.75})
          for (const double s : {0.3, 1.0 / 3.0, 0.75})
            {
              deallog << "Jacobi_homogenized_" << degree << "^(" << alpha << ','
                      << beta << ")(" << x << "," << s << ") ";

              // get the value
              const double value =
                jacobi_polynomial_homogenized_value(degree, alpha, beta, x, s);

              // compare to the definition Q_degree^(alpha,beta)(x,s) = s^degree
              // * P_degree^(alpha,beta)(x/s)
              const double definition =
                std::pow(s, degree) *
                jacobi_polynomial_value(degree, alpha, beta, x / s, true);

              if (std::abs(definition - value) < 1e-12)
                deallog << "ok ";
              else
                deallog << "value is: " << value
                        << " while the definition gives: " << definition
                        << " giving an error of  "
                        << std::abs(definition - value) << std::endl;

              deallog << std::endl;
            }
}



int
main()
{
  initlog();
  deallog.precision(10);

  test();
}
