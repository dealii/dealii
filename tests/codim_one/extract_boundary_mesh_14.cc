// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2018 by the deal.II authors
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

// check a bug in extract boundary mesh for a specific grid.

#include "../tests.h"

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

void
test()
{
  const int dim = 2;

  Triangulation<dim> triangulation;
  GridIn<dim>        gi;
  gi.attach_triangulation(triangulation);

  std::ifstream infile(SOURCE_DIR
                       "/../grid/grids/unstructured_refined_ball.msh");
  gi.read_msh(infile);

  // now extract the surface mesh
  Triangulation<dim - 1, dim> triangulation_surface;
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
