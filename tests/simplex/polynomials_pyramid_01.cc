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


// Like polyomials_01 for pyramids. Test PolynomialsPyramid on the points of the
// quadrature rule.


#include <deal.II/base/polynomials_pyramid.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_pyramid_p.h>

#include "../tests.h"


template <int dim>
void
test(const unsigned int degree)
{
  FE_PyramidP<dim> fe_pyramid(degree);
  const auto       support_points = fe_pyramid.get_unit_support_points();
  const auto       poly =
    ScalarLagrangePolynomialPyramid(degree,
                                    fe_pyramid.n_dofs_per_cell(),
                                    support_points);
  QGaussPyramid<dim> quad(2);

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
  {
    deallog.push("3d-3");
    test<3>(3);
    deallog.pop();
  }
}
