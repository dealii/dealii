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


// another failure that had to do that in the library we assumed that
// the left neighbor of the right neighbor of a cell is the cell
// itself. this holds true if dim==spacedim, but not
// otherwise. falsely making this assumption led to a strange failure
// in refine_global().

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

using namespace std;


void test ()
{
  const unsigned int spacedim = 2;
  const unsigned int dim = spacedim-1;

  Triangulation<dim,spacedim> boundary_mesh;
  map<Triangulation<dim,spacedim>::cell_iterator,Triangulation<spacedim,spacedim>::face_iterator >
  surface_to_volume_mapping;
  Triangulation<spacedim> volume_mesh;
  GridGenerator::hyper_cube(volume_mesh);
  volume_mesh.refine_global(1);
  surface_to_volume_mapping
    = GridGenerator::extract_boundary_mesh (volume_mesh, boundary_mesh);
  boundary_mesh.refine_global(1);

  for (Triangulation<dim,spacedim>::active_cell_iterator
       cell = boundary_mesh.begin_active();
       cell != boundary_mesh.end(); ++cell)
    {
      deallog << "Cell=" << cell << std::endl;
      deallog << "   neighbors: " << cell->neighbor(0) << ' ' << cell->neighbor(1)
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
