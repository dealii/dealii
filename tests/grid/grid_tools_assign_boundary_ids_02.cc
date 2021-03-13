// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Test assign_boundary_ids when overlapping predicates are provided.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"

void
test()
{
  Triangulation<3> tria;
  GridGenerator::hyper_cube(tria);
  // assign b-id of 3 to four faces where face_center[x] == 0.5
  // and b-id of 6 to two faces where
  //   face_center[x] + face_center[y] == 1.5
  // One of these faces already has b-id of 3, so we should not
  // change it. Therefore, only once face will contain b-id of 6.
  GridTools::assign_boundary_ids(
    tria,
    {{3, [](const auto &p) { return p[0] == 0.5; }},
     {6, [](const auto &p) { return p[0] + p[1] == 1.5; }}});
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      for (auto f : cell->face_indices())
        if (cell->face(f)->at_boundary())
          deallog << "face " << f << " B-id " << cell->face(f)->boundary_id()
                  << std::endl;
}

int
main()
{
  initlog();

  test();
  deallog << "OK" << std::endl;
  return 0;
}
