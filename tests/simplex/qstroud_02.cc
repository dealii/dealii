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

// Tests QStroudSimplex quadrature rules by comparison with QGaussSimplex.


#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"


template <int dim>
void
check_accuracy_2D(const unsigned int n_points_1D)
{
  const unsigned int accuracy = 2 * n_points_1D - 1;

  for (unsigned int i = 0; i < accuracy + 1; ++i)
    for (unsigned int j = 0; j + i < accuracy + 1; ++j)
      {
        Tensor<1, dim> monomial_powers;
        monomial_powers[0] = i;
        monomial_powers[1] = j;

        const Functions::Monomial<dim> func(monomial_powers);
        const QStroudSimplex<dim>      quad(n_points_1D);
        const QGaussSimplex<dim>       qgauss(n_points_1D);

        double integrand_stroud = 0.0;
        for (unsigned int q = 0; q < quad.size(); ++q)
          integrand_stroud += quad.weight(q) * func.value(quad.point(q));

        double integrand_gauss = 0.0;
        for (unsigned int q = 0; q < qgauss.size(); ++q)
          integrand_gauss += qgauss.weight(q) * func.value(qgauss.point(q));

        const double error = std::abs(integrand_stroud - integrand_gauss);

        deallog << "Monomial powers = " << monomial_powers << std::endl;

        if (error > 1e-14)
          {
            deallog << "With n points 1D " << n_points_1D << std::endl;
            deallog << "Integrand = " << integrand_stroud << std::endl;
            deallog << "Error = " << error << " Failed!" << std::endl;
          }
        else
          deallog << " passed!" << std::endl;
      }
}


template <int dim>
void
check_accuracy_3D(const unsigned int n_points_1D)
{
  const unsigned int accuracy = 2 * n_points_1D - 1;

  for (unsigned int i = 0; i < accuracy + 1; ++i)
    for (unsigned int j = 0; j + i < accuracy + 1; ++j)
      for (unsigned int k = 0; k + j + i < accuracy + 1; ++k)
        {
          Tensor<1, dim> monomial_powers;
          monomial_powers[0] = i;
          monomial_powers[1] = j;
          if constexpr (dim == 3)
            monomial_powers[2] = k;

          const Functions::Monomial<dim> func(monomial_powers);
          const QStroudSimplex<dim>      quad(n_points_1D);
          const QGaussSimplex<dim>       qgauss(n_points_1D);

          double integrand_stroud = 0.0;
          for (unsigned int q = 0; q < quad.size(); ++q)
            integrand_stroud += quad.weight(q) * func.value(quad.point(q));

          double integrand_gauss = 0.0;
          for (unsigned int q = 0; q < qgauss.size(); ++q)
            integrand_gauss += qgauss.weight(q) * func.value(qgauss.point(q));

          const double error = std::abs(integrand_stroud - integrand_gauss);

          deallog << "Monomial powers = " << monomial_powers << std::endl;

          if (error > 1e-14)
            {
              deallog << "With n points 1D " << n_points_1D << std::endl;
              deallog << "Integrand = " << integrand_stroud << std::endl;
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

  for (unsigned int i = 1; i < 6; ++i)
    {
      check_accuracy_2D<2>(i);
      check_accuracy_3D<3>(i);
    }
}
