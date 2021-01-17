// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// like _11, but refine the volume mesh and not the surface mesh. this
// should yield the same result

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold.h>
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

  for (const auto bid : triangulation.get_manifold_ids())
    if (bid != numbers::flat_manifold_id)
      triangulation.set_manifold(bid, FlatManifold<3>());

  static const CylindricalManifold<dim> outer_cylinder(0);
  triangulation.set_manifold(0, outer_cylinder);

  // refine the surface mesh to see the effect of boundary/manifold
  // indicators
  triangulation.refine_global(1);

  // now extract the surface mesh
  Triangulation<dim - 1, dim> triangulation_surface;

  for (const auto bid : triangulation.get_manifold_ids())
    if (bid != numbers::flat_manifold_id)
      triangulation_surface.set_manifold(bid, FlatManifold<2, 3>());

  static const CylindricalManifold<dim - 1, dim> surface_cyl(0);
  triangulation_surface.set_manifold(0, surface_cyl);

  GridGenerator::extract_boundary_mesh(triangulation, triangulation_surface);

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
