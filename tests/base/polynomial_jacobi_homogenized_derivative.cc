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

// check jacobi_polynomial_homogenized_derivative by comparing it to its
// definition

#include <deal.II/base/polynomial.h>

#include "../tests.h"

using namespace Polynomials;

void
test(const double tol)
{
  for (int alpha = 0; alpha < 3; ++alpha)
    for (int beta = 0; beta < 3; ++beta)
      for (unsigned int degree = 0; degree < 6; ++degree)
        for (const double x : {0.3, 1.0 / 3.0, 0.75})
          for (const double s : {0.3, 1.0 / 3.0, 0.75})
            {
              deallog << "Jacobi_homogenized_" << degree << "^(" << alpha << ','
                      << beta << ")(" << x << "," << s << ") ";

              // 0th derivative which should return the value
              const double deriv_0 = jacobi_polynomial_homogenized_derivative(
                0, 0, degree, alpha, beta, x, s);
              const double value =
                jacobi_polynomial_homogenized_value(degree, alpha, beta, x, s);

              if (std::abs(deriv_0 - value) < 1e-12)
                deallog << "ok ";
              else
                deallog << "0th derivative is: " << deriv_0
                        << " while value is: " << value
                        << " giving an error of  " << std::abs(deriv_0 - value)
                        << std::endl;


              // derivative in x
              const double deriv_x = jacobi_polynomial_homogenized_derivative(
                1, 0, degree, alpha, beta, x, s);
              // definition is Q_degree^(alpha,beta)(x,s) = s^degree *
              // P_degree^(alpha,beta)(x/s)
              const double deriv_definition_x =
                std::pow(s, degree) *
                jacobi_polynomial_derivative(degree, alpha, beta, x / s, true) /
                s;


              if (std::abs(deriv_x - deriv_definition_x) < 1e-12)
                deallog << "ok ";
              else
                deallog << "derivative in x is: " << deriv_x
                        << " while result from the definition is: "
                        << deriv_definition_x << " giving an error of  "
                        << std::abs(deriv_x - deriv_definition_x) << std::endl;

              // derivative in s
              const double deriv_s = jacobi_polynomial_homogenized_derivative(
                0, 1, degree, alpha, beta, x, s);
              const double deriv_definition_s =
                degree * std::pow(s, degree - 1) *
                  jacobi_polynomial_value(degree, alpha, beta, x / s, true) +
                std::pow(s, degree) *
                  jacobi_polynomial_derivative(
                    degree, alpha, beta, x / s, true) *
                  x / (s * s) * (-1.0);


              if (std::abs(deriv_s - deriv_definition_s) < 1e-12)
                deallog << "ok ";
              else
                deallog << "derivative in s is: " << deriv_s
                        << " while result from the definition is: "
                        << deriv_definition_s << " giving an error of  "
                        << std::abs(deriv_s - deriv_definition_s) << std::endl;

              deallog << std::endl;
            }
}



int
main()
{
  initlog();
  deallog.precision(10);

  const double tol = 1e-5;

  test(tol);
}
