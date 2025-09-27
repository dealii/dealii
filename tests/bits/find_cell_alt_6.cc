// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// On a 2D mesh of the following structure look for the cells surrounding
// each vertex, using the find_active_cell_around_point with Mapping:
//
// x-----x-----x
// |     |     |
// |     |     |
// |     |     |
// x--x--x-----x
// |  |  |     |
// x--x--x     x
// |  |  |     |
// x--x--x-----x

#include <deal.II/fe/mapping_q1.h>

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
  const std::vector<Point<2>> &v = tria.get_vertices();
  MappingQ<2>                  map(1);

  for (unsigned i = 0; i < tria.n_vertices(); ++i)
    {
      std::pair<Triangulation<2>::active_cell_iterator, Point<2>> cell =
        GridTools::find_active_cell_around_point(map, tria, v[i]);

      deallog << "Vertex <" << v[i] << "> found in cell ";
      for (const unsigned int v : GeometryInfo<2>::vertex_indices())
        deallog << '<' << cell.first->vertex(v) << "> ";
      deallog << " [local: " << cell.second << ']' << std::endl;
    }
}


int
main()
{
  initlog();

  try
    {
      Triangulation<2> coarse_grid;
      GridGenerator::hyper_cube(coarse_grid);
      coarse_grid.refine_global(1);
      coarse_grid.begin_active()->set_refine_flag();
      coarse_grid.execute_coarsening_and_refinement();
      check(coarse_grid);
    }
  catch (const std::exception &exc)
    {
      // we shouldn't get here...
      deallog << "Caught an error..." << std::endl;
      deallog << exc.what() << std::endl;
    }
}
