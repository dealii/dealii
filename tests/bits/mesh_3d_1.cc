//----------------------------  mesh_3d_1.cc  ---------------------------
//    mesh_3d_1.cc,v 1.2 2003/09/28 14:27:33 wolf Exp
//    Version: 
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_3d_1.cc  ---------------------------


// the real reason why coarsening_3d failed: take three cells to form an
// L, and we end up with an external face between two of them, which
// however has edges that are shared by the two cells. this creates a
// major upheaval in the data structures later on (as testified by
// coarsening_3d), but we can check this fact much earlier already (done
// here)

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
  std::ofstream logfile("mesh_3d_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  Triangulation<3> coarse_grid;
  create_L_shape (coarse_grid);

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
                << std::endl;
    }
}

  
  
