// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
