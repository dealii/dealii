// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Same as the first test, but on a 3D grid of the same structure

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
  for (unsigned i = 0; i < tria.n_vertices(); ++i)
    {
      std::vector<Triangulation<3>::active_cell_iterator> cells =
        GridTools::find_cells_adjacent_to_vertex(tria, i);

      deallog << "Vertex " << i << " at " << tria.get_vertices()[i] << ": "
              << cells.size() << " cells" << std::endl;

      for (unsigned c = 0; c < cells.size(); ++c)
        deallog << "   " << cells[c] << std::endl;
    }
}


int
main()
{
  initlog();

  try
    {
      Triangulation<3> coarse_grid;
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
