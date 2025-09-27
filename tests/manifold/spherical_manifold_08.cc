// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check SphericalManifold::get_intermediate_point for two particular points
// that are almost along a line to the center (this used to trigger when
// refining a 3d shell starting from 12 initial elements 6 times globally).

#include <deal.II/base/point.h>

#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(10);

  SphericalManifold<3> spherical;

  Point<3> p1(-0.036923837121444085,
              -0.87096345314118251,
              -0.075394254055706517);
  Point<3> p2(-0.037583190343442631,
              -0.8865163726762405,
              -0.076740572095662762);

  // when projected to radius one, both points are almost the same:
  // p1/p1.norm() : [-0.042198670995936091, -0.99538680358992271,
  // -0.086164861777950297] p2/p2.norm() : [-0.042198669859304011,
  // -0.99538680440841054, -0.086164852879340656]

  deallog << "intermediate point: "
          << spherical.get_intermediate_point(p1, p2, 0.5) << std::endl;

  return 0;
}
