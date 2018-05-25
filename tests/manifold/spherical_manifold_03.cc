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

// test tangent vectors to SphericalManifold at the poles

#include <deal.II/base/utilities.h>

#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"


int
main()
{
  initlog();

  const SphericalManifold<3> manifold;

  // get tangent vectors at the south pole of the sphere in direction
  // of the meridional equator point and a point 90 degrees to the
  // east or west of that point. this should yield two tangent vectors
  // that are orthogonal to each other
  deallog << manifold.get_tangent_vector(Point<3>(0, 0, -1), Point<3>(1, 0, 0))
          << std::endl
          << manifold.get_tangent_vector(Point<3>(0, 0, -1), Point<3>(0, 1, 0))
          << std::endl;
}
