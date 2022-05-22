// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"

// Test QIteratedSimplex accuracy

template <int dim>
void
print(const Quadrature<dim> &quad)
{
  deallog << "quad size = " << quad.size() << std::endl;
  for (unsigned int q = 0; q < quad.size(); ++q)
    {
      deallog << quad.point(q) << ' ';
      deallog << quad.weight(q) << ' ';
      deallog << std::endl;
    }
}

template <int dim>
void
check_accuracy_1D(const unsigned int n_points_1D, const unsigned int n_copies)
{
  const unsigned int accuracy = 2 * n_points_1D - 1;

  Tensor<1, dim> monomial_powers;
  unsigned int   sum = 0;
  for (unsigned int d = 0; d < dim; ++d)
    {
      monomial_powers[d] += accuracy / dim;
      sum += accuracy / dim;
    }

  // if we aren't at the correct degree then add the rest to the final
  // component
  monomial_powers[dim - 1] += accuracy - sum;

  const Functions::Monomial<dim> func(monomial_powers);
  const QIteratedSimplex<dim> quad(QWitherdenVincentSimplex<dim>(n_points_1D),
                                   n_copies);

  deallog << "Monomial powers = " << monomial_powers << std::endl;
  double integrand = 0.0;
  for (unsigned int q = 0; q < quad.size(); ++q)
    integrand += quad.weight(q) * func.value(quad.point(q));
  auto old_precision = deallog.precision(16);
  deallog << "Integrand = " << integrand << std::endl;
  deallog.precision(old_precision);
}

int
main()
{
  initlog();

  deallog << "2D:" << std::endl;
  print(QIteratedSimplex<2>(QWitherdenVincentSimplex<2>(2), 1));
  print(QIteratedSimplex<2>(QWitherdenVincentSimplex<2>(2), 2));
  deallog << std::endl << "3D:" << std::endl;
  print(QIteratedSimplex<3>(QWitherdenVincentSimplex<3>(2), 1));
  print(QIteratedSimplex<3>(QWitherdenVincentSimplex<3>(2), 2));

  deallog << std::endl << std::endl;
  check_accuracy_1D<2>(1, 1);
  check_accuracy_1D<2>(1, 2);
  check_accuracy_1D<2>(1, 4);
  check_accuracy_1D<2>(2, 1);
  check_accuracy_1D<2>(2, 2);
  check_accuracy_1D<2>(2, 4);
  check_accuracy_1D<2>(6, 1);
  check_accuracy_1D<2>(6, 2);
  check_accuracy_1D<2>(6, 4);

  check_accuracy_1D<3>(1, 1);
  check_accuracy_1D<3>(1, 2);
  check_accuracy_1D<3>(1, 4);
  check_accuracy_1D<3>(3, 1);
  check_accuracy_1D<3>(3, 2);
  check_accuracy_1D<3>(3, 4);
  check_accuracy_1D<3>(5, 1);
  check_accuracy_1D<3>(5, 2);
  check_accuracy_1D<3>(5, 4);
}
