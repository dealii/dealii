//----------------------------  mesh_3d_3.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_3d_3.cc  ---------------------------


// generate two cubes that are attached to each other in a way so that
// the edges are all ok, but the normals of the common face don't
// match up for the standard orientation of the normals. we thus have
// to store the face orientation in each cell
//
// for this grid, check that vertex numbers still match up

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
  std::ofstream logfile("mesh_3d_3.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Triangulation<3> coarse_grid;
  create_two_cubes (coarse_grid);

  const Triangulation<3>::active_cell_iterator
    cells[2] = { coarse_grid.begin_active(),
                 ++coarse_grid.begin_active() };

                                   // output all vertices
  for (unsigned int c=0; c<2; ++c)
    for (unsigned int v=0; v<GeometryInfo<3>::vertices_per_cell; ++v)
      deallog << "Cell " << c << ", vertex " << v
              << ": " << cells[c]->vertex_index(v)
              << "  @  " << cells[c]->vertex(v)
              << std::endl;

                                   // make sure by hand that certain
                                   // vertices match up
  Assert (cells[0]->vertex(1) == cells[1]->vertex(3),
          ExcInternalError());
  Assert (cells[0]->vertex(2) == cells[1]->vertex(2),
          ExcInternalError());
  Assert (cells[0]->vertex(5) == cells[1]->vertex(7),
          ExcInternalError());
  Assert (cells[0]->vertex(6) == cells[1]->vertex(6),
          ExcInternalError());
}

  
  
