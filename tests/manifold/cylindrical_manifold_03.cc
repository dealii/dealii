// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test cylindrical manifold for the mid point axis of a cylinder

#include "../tests.h"


// all include files you need here
#include <deal.II/grid/manifold_lib.h>



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  Tensor<1, 3>   axis({0.0, 0.0, 1.0});
  const Point<3> origin(1.0, 2.0, 3.0);

  // the larger the tolerance value, the more likely we identify the mid point
  const double tolerance = 1e-6;

  const CylindricalManifold<3> cylinder(axis, origin, tolerance);

  // take two points symmetric about cylinder axis
  const double          offset = 1.0;
  std::vector<Point<3>> surrounding_points_vector(
    {Point<3>(origin[0] + offset, origin[1], origin[2]),
     Point<3>(origin[0] - offset, origin[1], origin[2])});
  const ArrayView<const Point<3>> surrounding_points =
    make_array_view(surrounding_points_vector);
  const std::vector<double>     weights_vector({0.5, 0.5});
  const ArrayView<const double> weights = make_array_view(weights_vector);

  deallog << "New point is at "
          << cylinder.get_new_point(surrounding_points, weights)
          << " and it should be " << origin << std::endl;

  surrounding_points_vector[0][2] = -1;
  surrounding_points_vector[1][2] = -1;

  deallog << "New point is at "
          << cylinder.get_new_point(surrounding_points, weights)
          << " and it should be " << Point<3>(origin[0], origin[1], -1.0)
          << std::endl;

  axis[0] = 0.3;
  axis[1] = 0.6;
  axis[2] = 0.8;

  surrounding_points_vector[0][2] = origin[2];
  surrounding_points_vector[1][2] = origin[2];

  const CylindricalManifold<3> second_cylinder(axis, origin, tolerance);
  deallog << "New point is at "
          << second_cylinder.get_new_point(surrounding_points, weights)
          << " and it should be " << origin << std::endl;

  return 0;
}
