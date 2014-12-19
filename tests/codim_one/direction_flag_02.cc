// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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


// test direction flags in a 2d mesh embedded in 3d

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

using namespace std;


void test ()
{
  const unsigned int spacedim = 3;
  const unsigned int dim = spacedim-1;

  Triangulation<dim,spacedim> boundary_mesh;
  Triangulation<spacedim> volume_mesh;
  GridGenerator::hyper_cube(volume_mesh);
  GridGenerator::extract_boundary_mesh (volume_mesh, boundary_mesh);
  for (Triangulation<dim,spacedim>::active_cell_iterator
       cell = boundary_mesh.begin_active();
       cell != boundary_mesh.end(); ++cell)
    {
      deallog << "Cell=" << cell;
      deallog << ", direction flag="
              << (cell->direction_flag() ? "true" : "false")
              << std::endl;
    }

  boundary_mesh.refine_global(1);

  for (Triangulation<dim,spacedim>::active_cell_iterator
       cell = boundary_mesh.begin_active();
       cell != boundary_mesh.end(); ++cell)
    {
      deallog << "Cell=" << cell << std::endl;
      deallog << ", direction flag="
              << (cell->direction_flag() ? "true" : "false")
              << std::endl;
    }
}



int main ()
{
  ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test ();
}
