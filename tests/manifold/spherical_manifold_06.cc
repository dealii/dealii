// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check SphericalManifold for get_new_point issues at the origin.

#include <deal.II/base/utilities.h>

#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"

int
main()
{
  initlog();

  {
    Point<2>                      center(1.0, 0.0);
    const SphericalManifold<2, 2> manifold(center);
    Point<2>                      P1(1.0, 1.0);
    Point<2>                      P2(1.0, -1.0);
    std::vector<Point<2>>         points;
    std::vector<double>           weights;
    points.push_back(P1);
    points.push_back(P2);
    weights.push_back(0.5);
    weights.push_back(0.5);
    Point<2> mid_point = manifold.get_new_point(points, weights);
    deallog << mid_point << std::endl;
  }

  {
    Point<3>                      center(1.0, 0.0, 0.0);
    const SphericalManifold<3, 3> manifold(center);
    Point<3>                      P1(1.0, 1.0, 0.0);
    Point<3>                      P2(1.0, -1.0, 0.0);
    std::vector<Point<3>>         points;
    std::vector<double>           weights;
    points.push_back(P1);
    points.push_back(P2);
    weights.push_back(0.5);
    weights.push_back(0.5);
    Point<3> mid_point = manifold.get_new_point(points, weights);
    deallog << mid_point << std::endl;
  }

  return 0;
}
