// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that face orientation flags are properly inherited by
// counting them

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

#include "mesh_3d.h"



unsigned int
count_wrong_faces(const Triangulation<3> &tria)
{
  unsigned int count = 0;

  // count faces with "wrong"
  // orientation
  for (Triangulation<3>::active_cell_iterator cell = tria.begin_active();
       cell != tria.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<3>::face_indices())
      if (cell->face_orientation(f) == false)
        ++count;
  return count;
}



void
check(Triangulation<3> &tria)
{
  const unsigned int initial_count = count_wrong_faces(tria);
  for (unsigned int r = 0; r < 3; ++r)
    {
      tria.refine_global(1);
      const unsigned int count = count_wrong_faces(tria);
      deallog << "'Wrong' faces = " << count << std::endl;
      Assert(count == initial_count * (4 << (2 * r)), ExcInternalError());
    }

  {
    coarsen_global(tria);
    const unsigned int count = count_wrong_faces(tria);
    deallog << "'Wrong' faces = " << count << std::endl;
    Assert(count == initial_count * (4 << (2 * 1)), ExcInternalError());
  }

  {
    tria.refine_global(1);
    const unsigned int count = count_wrong_faces(tria);
    deallog << "'Wrong' faces = " << count << std::endl;
    Assert(count == initial_count * (4 << (2 * 2)), ExcInternalError());
  }
}



int
main()
{
  initlog();

  {
    Triangulation<3> coarse_grid;
    create_two_cubes(coarse_grid);
    check(coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    create_L_shape(coarse_grid);
    check(coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_ball(coarse_grid);
    check(coarse_grid);
  }
}
