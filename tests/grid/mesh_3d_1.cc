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



// the real reason why coarsening_3d failed: take three cells to form an
// L, and we end up with an external face between two of them, which
// however has edges that are shared by the two cells. this creates a
// major upheaval in the data structures later on (as testified by
// coarsening_3d), but we can check this fact much earlier already (done
// here)

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

#include "mesh_3d.h"



int
main()
{
  initlog();

  Triangulation<3> coarse_grid;
  create_L_shape(coarse_grid);

  // output all lines and faces
  for (Triangulation<3>::active_cell_iterator cell = coarse_grid.begin_active();
       cell != coarse_grid.end();
       ++cell)
    {
      deallog << "Cell = " << cell << std::endl;
      for (unsigned int i = 0; i < GeometryInfo<3>::lines_per_cell; ++i)
        deallog << "    Line = " << cell->line(i) << " : "
                << cell->line(i)->vertex_index(0) << " -> "
                << cell->line(i)->vertex_index(1) << std::endl;

      for (unsigned int i = 0; i < GeometryInfo<3>::quads_per_cell; ++i)
        deallog << "    Quad = " << cell->quad(i) << " : "
                << cell->quad(i)->vertex_index(0) << " -> "
                << cell->quad(i)->vertex_index(1) << " -> "
                << cell->quad(i)->vertex_index(2) << " -> "
                << cell->quad(i)->vertex_index(3) << std::endl;
    }
}
