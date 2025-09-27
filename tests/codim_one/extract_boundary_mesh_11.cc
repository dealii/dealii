// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// like _10, but this time try to use manifold ids that can be copied
// from the volume to the surface mesh

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



void
test()
{
  const int dim = 3;

  Triangulation<dim> triangulation;
  GridGenerator::cylinder(triangulation, 100, 200);

  // copy boundary indicators to manifold indicators for boundary
  // faces. for boundary zero (the outer hull of the cylinder), we
  // need to make sure that the adjacent edges are also all
  // correct. for the other boundaries, don't bother with adjacent
  // edges
  for (Triangulation<dim>::active_cell_iterator cell =
         triangulation.begin_active();
       cell != triangulation.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->face(f)->at_boundary())
        {
          if (cell->face(f)->boundary_id() == 0)
            cell->face(f)->set_all_manifold_ids(0);
          else
            cell->face(f)->set_manifold_id(cell->face(f)->boundary_id());
        }

  static const CylindricalManifold<dim> outer_cylinder(0);
  triangulation.set_manifold(0, outer_cylinder);

  // now extract the surface mesh
  Triangulation<dim - 1, dim> triangulation_surface;

  for (const auto bid : triangulation.get_manifold_ids())
    if (bid != numbers::flat_manifold_id)
      triangulation_surface.set_manifold(bid, FlatManifold<2, 3>());

  static const CylindricalManifold<dim - 1, dim> surface_cyl(0);
  triangulation_surface.set_manifold(0, surface_cyl);

  GridGenerator::extract_boundary_mesh(triangulation, triangulation_surface);

  // refine the surface mesh to see the effect of boundary/manifold
  // indicators
  triangulation_surface.refine_global(1);
  GridOut().write_gnuplot(triangulation_surface, deallog.get_file_stream());

  deallog << triangulation_surface.n_used_vertices() << std::endl;
  deallog << triangulation_surface.n_active_cells() << std::endl;
}


int
main()
{
  initlog();

  test();
}
