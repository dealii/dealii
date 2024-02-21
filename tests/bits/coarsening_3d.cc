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



// this test failed with an internal error somewhere in the coarsening
// functions

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

#include "../grid/mesh_3d.h"



int
main()
{
  initlog();

  Triangulation<3> coarse_grid;
  create_L_shape(coarse_grid);

  // refine once, then unrefine again
  coarse_grid.refine_global(1);
  for (Triangulation<3>::active_cell_iterator c = coarse_grid.begin_active();
       c != coarse_grid.end();
       ++c)
    c->set_coarsen_flag();
  coarse_grid.execute_coarsening_and_refinement();

  deallog << "ok." << std::endl;
}
