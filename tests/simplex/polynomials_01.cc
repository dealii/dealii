// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test BarycentricPolynomials on an the points of an arbitrary
// quadrature rule.


#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"


template <int dim>
void
test(const unsigned int degree)
{
  const auto    poly = BarycentricPolynomials<dim>::get_fe_p_basis(degree);
  QSimplex<dim> quad(QGauss<dim>(degree + 1));

  std::vector<double>         values(poly.n());
  std::vector<Tensor<1, dim>> grads(poly.n());
  std::vector<Tensor<2, dim>> grad_grads;
  std::vector<Tensor<3, dim>> third_derivatives;
  std::vector<Tensor<4, dim>> fourth_derivatives;

  for (unsigned int i = 0; i < quad.size(); ++i)
    {
      poly.evaluate(quad.point(i),
                    values,
                    grads,
                    grad_grads,
                    third_derivatives,
                    fourth_derivatives);

      for (auto v : values)
        deallog << v << ' ';
      deallog << std::endl;

      for (auto v : grads)
        deallog << v << ' ';
      deallog << std::endl;
    }
}

int
main()
{
  initlog();

  {
    deallog.push("1d-1");
    test<1>(1);
    deallog.pop();
  }
  {
    deallog.push("1d-2");
    test<1>(2);
    deallog.pop();
  }
  {
    deallog.push("2d-1");
    test<2>(1);
    deallog.pop();
  }
  {
    deallog.push("2d-2");
    test<2>(2);
    deallog.pop();
  }
  {
    deallog.push("3d-1");
    test<3>(1);
    deallog.pop();
  }
  {
    deallog.push("3d-2");
    test<3>(2);
    deallog.pop();
  }
}
