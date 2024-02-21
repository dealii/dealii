// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// like find_cell_2, but with the strange meshes from the mesh_3d_* tests

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

#include "../grid/mesh_3d.h"



void
check(Triangulation<3> &tria)
{
  Point<3> p(1. / 3., 1. / 2., -1. / 5.);

  Triangulation<3>::active_cell_iterator cell =
    GridTools::find_active_cell_around_point(tria, p);

  deallog << cell << std::endl;
  for (const unsigned int v : GeometryInfo<3>::vertex_indices())
    deallog << '<' << cell->vertex(v) << "> ";
  deallog << std::endl;

  Assert(p.distance(cell->center()) < cell->diameter() / 2, ExcInternalError());
}


int
main()
{
  initlog();

  {
    Triangulation<3> coarse_grid;
    create_two_cubes(coarse_grid);
    coarse_grid.refine_global(1);
    check(coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    create_L_shape(coarse_grid);
    coarse_grid.refine_global(1);
    check(coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_ball(coarse_grid);
    coarse_grid.reset_all_manifolds();
    coarse_grid.refine_global(1);
    check(coarse_grid);
  }
}
