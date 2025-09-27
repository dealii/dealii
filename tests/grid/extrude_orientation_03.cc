// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test GridGenerator::extrude_triangulation taking slice z-coordinate values.
// This test is just a replicate of the one in extrude_orientation_02.cc except
// that the newly created overload of GridGenerator::extrude_triangulation is
// tested for face and edge orientations for the resulting 3d mesh.
// See https://github.com/dealii/dealii/issues/6158.
//
// test this for a circle extruded to a cylinder

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

void
test()
{
  Triangulation<2> tr;
  GridGenerator::hyper_ball(tr);

  for (Triangulation<2>::active_cell_iterator c = tr.begin_active();
       c != tr.end();
       ++c)
    {
      deallog << "2d cell " << c
              << " has the following face orientations:" << std::endl;
      for (const unsigned int l : GeometryInfo<2>::face_indices())
        deallog << "    " << (c->face_orientation(l) ? "true" : "false")
                << std::endl;
    }

  Triangulation<3>    tr3;
  std::vector<double> slice_points = {0, 0.1, 0.5};
  GridGenerator::extrude_triangulation(tr, slice_points, tr3);

  for (Triangulation<3>::active_cell_iterator c = tr3.begin_active();
       c != tr3.end();
       ++c)
    {
      deallog
        << "3d cell " << c
        << " has the following face orientation/flips and edge orientations:"
        << std::endl;
      for (const unsigned int f : GeometryInfo<3>::face_indices())
        deallog << "    face=" << f
                << (c->face_orientation(f) ? " -> true" : " -> false")
                << (c->face_flip(f) ? "/true" : "/false") << std::endl;
      for (unsigned int e = 0; e < GeometryInfo<3>::lines_per_cell; ++e)
        deallog << "    edge=" << e
                << (c->line_orientation(e) ==
                        numbers::default_geometric_orientation ?
                      " -> true" :
                      " -> false")
                << std::endl;
    }
}


int
main()
{
  initlog();

  test();
}
