// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Tests QGaussPyramid by integrating monomials over the reference element
// and comparing the result with the analytical solution.

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"

unsigned int
factorial(const unsigned int n)
{
  if (n < 2)
    return 1;
  return factorial(n - 1) * n;
}


template <int dim>
void
check_accuracy_1D(const unsigned int n_points_1D)
{
  const unsigned int accuracy = 2 * n_points_1D - 1;

  // Pyramid has Gauss quadrature in x,y and Gauss-Jacobi in z direction
  // mononomial leads to z^k*(1-z)^(i+j+2)

  for (const unsigned int i :
       std::vector<unsigned int>{{0, 1, 2, accuracy + 1}})
    for (const unsigned int j :
         std::vector<unsigned int>{{0, 1, 2, accuracy + 1}})
      for (unsigned int k = 0; i + j + k < accuracy - 1; ++k)
        {
          Tensor<1, dim> monomial_powers;
          monomial_powers[0] = i;
          monomial_powers[1] = j;
          monomial_powers[2] = k;

          const Functions::Monomial<dim> func(monomial_powers);
          const QGaussPyramid<dim>       quad(n_points_1D);

          double integrand = 0.0;
          for (unsigned int q = 0; q < quad.size(); ++q)
            integrand += quad.weight(q) * func.value(quad.point(q));

          const double analytical_solution =
            (i % 2 == 1 || j % 2 == 1) ?
              0.0 :
              4.0 / ((i + 1) * (j + 1)) * factorial(k) * factorial(i + j + 2) /
                factorial(i + j + k + 3);

          const double error = std::abs(integrand - analytical_solution);

          deallog << "Monomial powers = " << monomial_powers << std::endl;

          if (error > 1e-14)
            {
              deallog << "With n points 1D " << n_points_1D << std::endl;
              deallog << "Integrand = " << integrand << std::endl;
              deallog << "Analytical_solution = " << analytical_solution
                      << std::endl;
              deallog << "Error = " << error << " Failed!" << std::endl;
            }
          else
            deallog << " passed!" << std::endl;
        }
}



int
main()
{
  initlog();

  for (unsigned int i = 1; i < 7; ++i)
    check_accuracy_1D<3>(i);
}
