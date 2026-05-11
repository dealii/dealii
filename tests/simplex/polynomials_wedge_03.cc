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


// Test PolynomialsWedge by checking if the basis sums up to 1 at arbitrary
// points and if the gradient is 0 at any point


#include <deal.II/base/polynomials_wedge.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_wedge_p.h>

#include "../tests.h"


template <int dim>
void
test(const unsigned int degree)
{
  FE_WedgeP<dim> fe_wedge(degree);
  auto           support_points = fe_wedge.get_unit_support_points();
  const auto     poly = ScalarLagrangePolynomialWedge(degree, support_points);

  double         sum;
  Tensor<1, dim> grad_sum;

  std::vector<Point<dim>> points;
  for (const auto &p : support_points)
    points.emplace_back(p);

  // Add some random points
  points.emplace_back(Point<dim>(0.36658, 0.58775, 0.21455));
  points.emplace_back(Point<dim>(0.64464, 0.3546, 0.3246));

  // Add points from the quadrature
  QGaussWedge<dim> quad(2);
  for (unsigned int i = 0; i < quad.size(); ++i)
    points.emplace_back(quad.point(i));

  for (const auto &p : points)
    {
      sum      = 0.;
      grad_sum = 0.;
      for (unsigned int j = 0; j < support_points.size(); ++j)
        {
          sum += poly.compute_value(j, p);
          grad_sum += poly.compute_grad(j, p);
        }
      if (std::abs(sum - 1) < 1e-12)
        deallog << "ok" << std::endl;
      else
        deallog << "Point " << p << " not ok " << sum << std::endl;

      for (unsigned int d = 0; d < dim; ++d)
        {
          if (std::abs(grad_sum[d]) < 1e-12)
            deallog << "ok ";
          else
            deallog << "Point " << p << " gradient not ok " << grad_sum[d];
        }
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
    deallog.push("3d-2");
    test<3>(2);
    deallog.pop();
  }
}
