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

// check jacobi_polynomial_derivative by approximating the derivative with a
// finite difference stencil

#include <deal.II/base/polynomial.h>

#include "../tests.h"

using namespace Polynomials;


double
get_finite_difference_first_derivative(const double       h,
                                       const unsigned int degree,
                                       const int          alpha,
                                       const int          beta,
                                       const double       x,
                                       const bool         rescale)
{
  const double deriv =
    (jacobi_polynomial_value(degree, alpha, beta, x + h, rescale) -
     jacobi_polynomial_value(degree, alpha, beta, x - h, rescale)) /
    (2.0 * h);
  return deriv;
}

double
get_finite_difference_second_derivative(const double       h,
                                        const unsigned int degree,
                                        const int          alpha,
                                        const int          beta,
                                        const double       x,
                                        const bool         rescale)
{
  const double deriv =
    (jacobi_polynomial_value(degree, alpha, beta, x + h, rescale) -
     2.0 * jacobi_polynomial_value(degree, alpha, beta, x, rescale) +
     jacobi_polynomial_value(degree, alpha, beta, x - h, rescale)) /
    (h * h);
  return deriv;
}

void
test(const double tol)
{
  for (int alpha = 0; alpha < 3; ++alpha)
    for (int beta = 0; beta < 3; ++beta)
      for (unsigned int degree = 0; degree < 6; ++degree)
        for (const double x : {0.3, 1.0 / 3.0, 0.75})
          for (const bool rescale : {true, false})
            {
              deallog << "Jacobi_" << degree << "^(" << alpha << ',' << beta
                      << ")(" << x << ") ";

              // 0th derivative which should return the value
              const double deriv_0 = jacobi_polynomial_kth_derivative(
                0, degree, alpha, beta, x, rescale);
              const double value =
                jacobi_polynomial_value(degree, alpha, beta, x, rescale);

              if (std::abs(deriv_0 - value) < tol)
                deallog << "ok ";
              else
                deallog << "0th derivative is: " << deriv_0
                        << " while value is: " << value
                        << " giving an error of  " << std::abs(deriv_0 - value)
                        << std::endl;


              // first derivative
              const double deriv =
                jacobi_polynomial_derivative(degree, alpha, beta, x, rescale);
              const double deriv_approx =
                get_finite_difference_first_derivative(
                  tol, degree, alpha, beta, x, rescale);

              if (std::abs(deriv - deriv_approx) < tol)
                deallog << "ok ";
              else
                deallog << "derivative is: " << deriv
                        << " while finite difference is: " << deriv_approx
                        << " giving an error of  "
                        << std::abs(deriv - deriv_approx) << std::endl;

              // second derivative
              const double second_deriv = jacobi_polynomial_kth_derivative(
                2, degree, alpha, beta, x, rescale);
              const double second_deriv_approx =
                get_finite_difference_second_derivative(
                  tol, degree, alpha, beta, x, rescale);

              if (std::abs(second_deriv - second_deriv_approx) < tol)
                deallog << "ok ";
              else
                deallog << "2nd derivative is: " << second_deriv
                        << " while finite difference is: "
                        << second_deriv_approx << " giving an error of "
                        << std::abs(second_deriv - second_deriv_approx)
                        << std::endl;
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
