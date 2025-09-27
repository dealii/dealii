// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test GridGenerator::replicate_triangulation's ability to copy boundary ids
// correctly

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

void
test()
{
  // no subcell data in 1D so nothing to test

  {
    Triangulation<2, 2> tr2;
    GridGenerator::hyper_cube(tr2, 0.0, 1.0, true);
    const std::vector<unsigned int> reps2 = {{2, 2}};
    Triangulation<2, 2>             res2;
    GridGenerator::replicate_triangulation(tr2, reps2, res2);
    for (const auto &face : res2.active_face_iterators())
      if (face->at_boundary())
        deallog << "face: " << face->center()
                << " boundary id: " << face->boundary_id() << std::endl;
  }

  {
    Triangulation<3, 3> tr3;
    GridGenerator::hyper_cube(tr3, 0.0, 1.0, true);
    // actually set some boundary line boundary ids to something besides zero
    // for testing purposes
    for (auto &face : tr3.active_face_iterators())
      {
        if (face->boundary_id() == 4)
          face->set_all_boundary_ids(4);
        else if (face->boundary_id() == 5)
          face->set_all_boundary_ids(5);
      }
    // set all (i.e., overwrite two set above) lines on the x = 1 face to
    // boundary id 1
    for (auto &face : tr3.active_face_iterators())
      if (face->boundary_id() == 1)
        face->set_all_boundary_ids(1);
    std::vector<unsigned int> reps3 = {2, 1, 1};
    Triangulation<3, 3>       res3;
    GridGenerator::replicate_triangulation(tr3, reps3, res3);
    for (const auto &face : res3.active_face_iterators())
      if (face->at_boundary())
        {
          deallog << "face: " << face->center()
                  << " boundary id: " << face->boundary_id() << std::endl;
          for (unsigned int l = 0; l < GeometryInfo<2>::lines_per_cell; ++l)
            deallog << "line: " << face->line(l)->center()
                    << " boundary id: " << face->line(l)->boundary_id()
                    << std::endl;
        }
  }
}

int
main()
{
  initlog();
  test();
}
