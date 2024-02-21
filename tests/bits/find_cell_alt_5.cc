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



// take a 3d mesh and check that we can find an arbitrary point's cell
// in it.

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


void
check(Triangulation<3> &tria)
{
  MappingQ<3> map(1);
  Point<3>    p(0.75, 0.75, 0.75);

  std::pair<Triangulation<3>::active_cell_iterator, Point<3>> cell =
    GridTools::find_active_cell_around_point(map, tria, p);

  deallog << cell.first << std::endl;
  for (const unsigned int v : GeometryInfo<3>::vertex_indices())
    deallog << '<' << cell.first->vertex(v) << "> ";
  deallog << "[ " << cell.second << "] ";
  deallog << std::endl;

  Assert(GeometryInfo<3>::distance_to_unit_cell(cell.second) < 1e-10,
         ExcInternalError());
}


int
main()
{
  initlog();

  try
    {
      Triangulation<3> coarse_grid;
      GridGenerator::hyper_cube(coarse_grid);
      coarse_grid.refine_global(3);
      check(coarse_grid);
    }
  catch (const std::exception &exc)
    {
      // we shouldn't get here...
      deallog << "Caught an error..." << std::endl;
      deallog << exc.what() << std::endl;
    }
}
