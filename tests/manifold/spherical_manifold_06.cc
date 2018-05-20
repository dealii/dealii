// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Check SphericalManifold for get_new_point issues at the origin.

#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/grid/manifold_lib.h>

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
