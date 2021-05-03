// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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



// check that face orientation flags are properly inherited by
// looking at them and checking their children

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

#include "mesh_3d.h"



void check_this(Triangulation<3> &tria)
{
  // look at all faces, not only
  // active ones
  for (Triangulation<3>::cell_iterator cell = tria.begin(); cell != tria.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<3>::face_indices())
      if (cell->has_children())
        for (unsigned int c = 0; c < GeometryInfo<3>::max_children_per_face;
             ++c)
          {
            Assert(cell->face_orientation(f) ==
                     cell
                       ->child(GeometryInfo<3>::child_cell_on_face(
                         RefinementCase<3>::isotropic_refinement, f, c))
                       ->face_orientation(f),
                   ExcInternalError());
            deallog << "Cell << " << cell << ", face " << f << " subface " << c
                    << " is ok." << std::endl;
          }
}


void check(Triangulation<3> &tria)
{
  deallog << "Initial check" << std::endl;
  check_this(tria);

  for (unsigned int r = 0; r < 3; ++r)
    {
      tria.refine_global(1);
      deallog << "Check " << r << std::endl;
      check_this(tria);
    }

  coarsen_global(tria);
  deallog << "Check " << 1 << std::endl;
  check_this(tria);

  tria.refine_global(1);
  deallog << "Check " << 2 << std::endl;
  check_this(tria);
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
