//----------------------------  mesh_3d_2.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_3d_2.cc  ---------------------------


// generate two cubes that are attached to each other in a way so that
// the edges are all ok, but the normals of the common face don't
// match up for the standard orientation of the normals. we thus have
// to store the face orientation in each cell
//
// this test just checks that such a mesh can be created. further
// tests check that we are still producing consistent states with this

#include "../tests.h"
#include "mesh_3d.h"

#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_reordering.h>

#include <fstream>


int main () 
{
  std::ofstream logfile("mesh_3d_2.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

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
  Assert ((++coarse_grid.begin_active())->face_orientation(4)
          == false,
          ExcInternalError());
}

  
  
