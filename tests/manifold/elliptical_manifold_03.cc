// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check gradient of elliptical manifold

#include "../tests.h"


// all include files you need here
#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <vector>



int
main()
{
  initlog();

  deallog << std::setprecision(10);

  // setup elliptical manifold
  const Tensor<1, 2> major_axis({2.0, 1.0});
  const Point<2>     center(7.0, 6.0);
  const double       eccentricity = 0.1;

  const EllipticalManifold<2> ellipse(center, major_axis, eccentricity);

  // reference point for computing gradient
  const Point<2> chart_point(1.0, 0.25 * numbers::PI);

  // compute gradient by central finite differences
  const double   h = 1e-6;
  const Point<2> chart_point_plusc(chart_point[0] + h, chart_point[1]);
  const Point<2> chart_point_minusc(chart_point[0] - h, chart_point[1]);
  const Point<2> chart_point_plusphi(chart_point[0], chart_point[1] + h);
  const Point<2> chart_point_minusphi(chart_point[0], chart_point[1] - h);

  const Point<2> space_point_plusc  = ellipse.push_forward(chart_point_plusc);
  const Point<2> space_point_minusc = ellipse.push_forward(chart_point_minusc);
  const Point<2> space_point_plusphi =
    ellipse.push_forward(chart_point_plusphi);
  const Point<2> space_point_minusphi =
    ellipse.push_forward(chart_point_minusphi);

  deallog << "Gradient by finite differences: "
          << (space_point_plusc - space_point_minusc) / (2.0 * h) << ' '
          << (space_point_plusphi - space_point_minusphi) / (2.0 * h)
          << std::endl;
  deallog << "Analytic gradient:              "
          << transpose(Tensor<2, 2>(ellipse.push_forward_gradient(chart_point)))
          << std::endl;

  return 0;
}
