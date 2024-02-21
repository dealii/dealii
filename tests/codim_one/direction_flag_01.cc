// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test direction flags in a 1d mesh embedded in 2d

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


void
test()
{
  const unsigned int spacedim = 2;
  const unsigned int dim      = spacedim - 1;

  Triangulation<dim, spacedim> boundary_mesh;
  Triangulation<spacedim>      volume_mesh;
  GridGenerator::hyper_cube(volume_mesh);
  GridGenerator::extract_boundary_mesh(volume_mesh, boundary_mesh);
  for (Triangulation<dim, spacedim>::active_cell_iterator cell =
         boundary_mesh.begin_active();
       cell != boundary_mesh.end();
       ++cell)
    {
      deallog << "Cell=" << cell;
      deallog << ", direction flag="
              << (cell->direction_flag() ? "true" : "false") << std::endl;
    }

  boundary_mesh.refine_global(1);

  for (Triangulation<dim, spacedim>::active_cell_iterator cell =
         boundary_mesh.begin_active();
       cell != boundary_mesh.end();
       ++cell)
    {
      deallog << "Cell=" << cell << std::endl;
      deallog << ", direction flag="
              << (cell->direction_flag() ? "true" : "false") << std::endl;
    }
}



int
main()
{
  initlog();

  test();
}
