// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// another failure that had to do that in the library we assumed that
// the left neighbor of the right neighbor of a cell is the cell
// itself. this holds true if dim==spacedim, but not
// otherwise. falsely making this assumption led to a strange failure
// in refine_global().

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
  std::map<Triangulation<dim, spacedim>::cell_iterator,
           Triangulation<spacedim, spacedim>::face_iterator>
                          surface_to_volume_mapping;
  Triangulation<spacedim> volume_mesh;
  GridGenerator::hyper_cube(volume_mesh);
  volume_mesh.refine_global(1);
  surface_to_volume_mapping =
    GridGenerator::extract_boundary_mesh(volume_mesh, boundary_mesh);
  boundary_mesh.refine_global(1);

  for (Triangulation<dim, spacedim>::active_cell_iterator cell =
         boundary_mesh.begin_active();
       cell != boundary_mesh.end();
       ++cell)
    {
      deallog << "Cell=" << cell << std::endl;
      deallog << "   neighbors: " << cell->neighbor(0) << ' '
              << cell->neighbor(1) << std::endl;
    }
}



int
main()
{
  initlog();

  test();
}
