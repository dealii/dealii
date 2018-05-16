// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// like _11, but using a different way of copying manifold ids

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>



void
test()
{
  const int dim=3;

  Triangulation<dim>   triangulation;
  GridGenerator::cylinder(triangulation, 100, 200);

  // copy boundary indicators to manifold indicators for boundary
  // faces
  GridTools::copy_boundary_to_manifold_id (triangulation);

  static const CylindricalManifold<dim> outer_cylinder (0);
  triangulation.set_manifold(0,outer_cylinder);

  // now extract the surface mesh
  Triangulation<dim-1,dim> triangulation_surface;

  static const CylindricalManifold<dim-1,dim> surface_cyl(0);
  triangulation_surface.set_manifold(0,surface_cyl);

  GridGenerator::extract_boundary_mesh(triangulation,triangulation_surface);

  // refine the surface mesh to see the effect of boundary/manifold
  // indicators
  triangulation_surface.refine_global (1);
  GridOut().write_gnuplot(triangulation_surface, deallog.get_file_stream());

  deallog << triangulation_surface.n_used_vertices() << std::endl;
  deallog << triangulation_surface.n_active_cells() << std::endl;
}


int
main ()
{
  initlog();
  deallog.depth_console(0);

  test();
}
