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


// Like polyomials_01 for simplices. Test PolynomialsSimplex on the points of
// the quadrature rule.


#include <deal.II/base/polynomials_simplex.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_simplex_p.h>

#include "../tests.h"


template <int dim>
void
test(const unsigned int degree)
{
  FE_SimplexP<dim> fe_simplex(degree);
  const auto       support_points = fe_simplex.get_unit_support_points();
  const auto poly = ScalarLagrangePolynomialSimplex(degree, support_points);
  QGaussSimplex<dim> quad(2);

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

  for (unsigned int i = 1; i < 4; ++i)
    {
      {
        deallog.push("1d-" + std::to_string(i));
        test<1>(i);
        deallog.pop();
      }
      {
        deallog.push("2d-" + std::to_string(i));
        test<2>(i);
        deallog.pop();
      }
      {
        deallog.push("3d-" + std::to_string(i));
        test<3>(i);
        deallog.pop();
      }
    }
}
