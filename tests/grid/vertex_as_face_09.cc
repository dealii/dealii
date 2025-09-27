// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// verify that we can do things like cell->face() in 1d as well. here:
// Triangulation::get_boundary_ids() should return an empty vector when
// called to create a mesh that is a loop


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



void
test()
{
  Triangulation<2, 2> volume_mesh;
  GridGenerator::hyper_cube(volume_mesh);

  Triangulation<1, 2> tria;

  GridGenerator::extract_boundary_mesh(volume_mesh, tria);

  deallog << "n_cells = " << tria.n_active_cells() << std::endl;
  deallog << "n_boundary_ids = " << tria.get_boundary_ids().size() << std::endl;
}



int
main()
{
  initlog();

  test();

  return 0;
}
