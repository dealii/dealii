// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// generate two cubes that are attached to each other in a way so that
// the edges are all ok, but the normals of the common face don't
// match up for the standard orientation of the normals. we thus have
// to store the face orientation in each cell
//
// this test just checks that such a mesh can be created. further
// tests check that we are still producing consistent states with this

#include "../tests.h"
#include "mesh_3d.h"

#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_reordering.h>

#include <fstream>


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Triangulation<3> coarse_grid;
  create_two_cubes (coarse_grid);

  // output all lines and faces
  for (Triangulation<3>::active_cell_iterator cell=coarse_grid.begin_active();
       cell != coarse_grid.end(); ++cell)
    {
      deallog << "Cell = " << cell << std::endl;
      for (unsigned int i=0; i<GeometryInfo<3>::lines_per_cell; ++i)
        deallog << "    Line = " << cell->line(i)
                << " : " << cell->line(i)->vertex_index(0)
                << " -> " << cell->line(i)->vertex_index(1)
                << std::endl;

      for (unsigned int i=0; i<GeometryInfo<3>::quads_per_cell; ++i)
        deallog << "    Quad = " << cell->quad(i)
                << " : " << cell->quad(i)->vertex_index(0)
                << " -> " << cell->quad(i)->vertex_index(1)
                << " -> " << cell->quad(i)->vertex_index(2)
                << " -> " << cell->quad(i)->vertex_index(3)
                << std::endl
                << "           orientation = "
                << (cell->face_orientation(i) ? "true" : "false")
                << std::endl;
    }

  // we know that from the second
  // cell, the common face must have
  // wrong orientation. check this
  Assert ((++coarse_grid.begin_active())->face_orientation(5)
          == false,
          ExcInternalError());
}



