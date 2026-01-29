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

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"

template <int dim>
void
test(const unsigned int n_points_1D)
{
  deallog << "N points 1D: " << n_points_1D << std::endl;
  std::vector<Point<dim>> quadrature_points;
  std::vector<double>     weights;

  if (n_points_1D == 1)
    {
      const double Q14 = 1.0 / 4.0;
      const double Q43 = 4.0 / 3.0;

      quadrature_points.emplace_back(0, 0, Q14);
      weights.emplace_back(Q43);
    }
  else if (n_points_1D == 2)
    {
      quadrature_points.emplace_back(-0.50661630334979,
                                     -0.50661630334979,
                                     0.12251482265544);
      quadrature_points.emplace_back(-0.26318405556971,
                                     -0.26318405556971,
                                     0.54415184401122);
      quadrature_points.emplace_back(-0.50661630334979,
                                     +0.50661630334979,
                                     0.12251482265544);
      quadrature_points.emplace_back(-0.26318405556971,
                                     +0.26318405556971,
                                     0.54415184401122);
      quadrature_points.emplace_back(+0.50661630334979,
                                     -0.50661630334979,
                                     0.12251482265544);
      quadrature_points.emplace_back(+0.26318405556971,
                                     -0.26318405556971,
                                     0.54415184401122);
      quadrature_points.emplace_back(+0.50661630334979,
                                     +0.50661630334979,
                                     0.12251482265544);
      quadrature_points.emplace_back(+0.26318405556971,
                                     +0.26318405556971,
                                     0.54415184401122);

      weights.emplace_back(0.23254745125351);
      weights.emplace_back(0.10078588207983);
      weights.emplace_back(0.23254745125351);
      weights.emplace_back(0.10078588207983);
      weights.emplace_back(0.23254745125351);
      weights.emplace_back(0.10078588207983);
      weights.emplace_back(0.23254745125351);
      weights.emplace_back(0.10078588207983);
    }

  QGaussPyramid<dim> quad(n_points_1D);
  deallog << "quad size " << quad.size() << " "
          << " reference size " << weights.size() << std::endl;

  for (unsigned int i = 0; i < quad.size(); ++i)
    {
      deallog << quad.point(i) << " " << quad.weight(i) << std::endl;
      const auto p = quad.point(i);
      if (std::abs(p.distance_square(quadrature_points[i])) < 1e-12)
        deallog << "Correct point" << std::endl;
      else
        deallog << "Point does not match" << std::endl;

      if (std::abs(quad.weight(i) - weights[i]) < 1e-12)
        deallog << "Correct weight" << std::endl;
      else
        deallog << "Weight does not match" << std::endl;
    }

  deallog << std::endl;
}


int
main()
{
  initlog();

  test<3>(1);
  test<3>(2);
}
