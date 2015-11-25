// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2015 by the deal.II authors
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


// We can not copy boundary ids from edges of the volume mesh to edges
// of the surface mesh in 3d because the boundary edges of the volume
// mesh are not necessarily boundary edges of the surface mesh.
//
// if the original mesh's geometry was only defined via the boundary
// indicators, this leads to a surface mesh with a geometry
// description of edges that is not likely what was desired. there is
// little one can do about it with only boundary_ids, and without
// resorting to manifold_ids. this testcase in essence just ensures
// that we continue to get this "undesirable" behavior until someone
// comes up with a better idea

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_boundary_lib.h>



void test()
{
  const int dim=3;

  Triangulation<dim>   triangulation;
  GridGenerator::cylinder(triangulation, 100, 200);

  static const CylinderBoundary<dim> outer_cylinder (100,0);
  triangulation.set_boundary(0,outer_cylinder);

  // now extract the surface mesh
  Triangulation<dim-1,dim> triangulation_surface;

  static const CylinderBoundary<dim-1,dim> surface_cyl(100,0);
  triangulation_surface.set_boundary(0,surface_cyl);

  GridGenerator::extract_boundary_mesh(triangulation,triangulation_surface);

  // refine the surface mesh to see the effect of boundary/manifold
  // indicators
  triangulation_surface.refine_global (1);
  GridOut().write_gnuplot(triangulation_surface, deallog.get_file_stream());

  deallog << triangulation_surface.n_used_vertices() << std::endl;
  deallog << triangulation_surface.n_active_cells() << std::endl;
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test();
}
