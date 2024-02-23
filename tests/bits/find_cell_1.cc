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



// take a 2d mesh and check that we can find an arbitrary point's cell
// in it

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



void
check(Triangulation<2> &tria)
{
  Point<2> p(1. / 3., 1. / 2. - 1e-10); // avoid ambiguity for hypercube mesh

  Triangulation<2>::active_cell_iterator cell =
    GridTools::find_active_cell_around_point(tria, p);

  deallog << cell << std::endl;
  for (const unsigned int v : GeometryInfo<2>::vertex_indices())
    deallog << '<' << cell->vertex(v) << "> ";
  deallog << std::endl;

  Assert(p.distance(cell->center()) < cell->diameter() / 2, ExcInternalError());
}


int
main()
{
  initlog();

  {
    Triangulation<2> coarse_grid;
    GridGenerator::hyper_cube(coarse_grid);
    coarse_grid.refine_global(2);
    check(coarse_grid);
  }

  {
    Triangulation<2> coarse_grid;
    GridGenerator::hyper_ball(coarse_grid);
    static const SphericalManifold<2> boundary;
    coarse_grid.set_manifold(0, boundary);
    coarse_grid.refine_global(2);
    check(coarse_grid);
  }
}
