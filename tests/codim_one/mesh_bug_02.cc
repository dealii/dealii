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


// Test that material_id and manifold_id don't interfere with each other.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"


int
main()
{
  initlog();

  // triangulation is a single 2d square element embedded into 3d space
  Triangulation<2, 3>   plane;
  Point<3>              p1(1., 0., 0.);
  Point<3>              p2(1., 1., 0.);
  Point<3>              p3(1., 0., 1.);
  Point<3>              p4(1., 1., 1.);
  std::vector<Point<3>> vertices = {p1, p2, p3, p4};
  GridGenerator::general_cell(plane, vertices);

  // set the material id of the cell to 1
  plane.begin_active()->set_material_id(1);

  // now associate a spherical manifold object with manifold_id 1
  //(this object isn't used in this minimal example, but there could
  // in principle be other cells to which the manifold_id 1 is assigned)
  SphericalManifold<2, 3> spherical_manifold;
  plane.set_manifold(1, spherical_manifold);

  plane.refine_global(2);

  GridOut gridout;
  gridout.write_msh(plane, deallog.get_file_stream());
}
