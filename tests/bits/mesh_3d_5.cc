//----------------------------  mesh_3d_5.cc  ---------------------------
//    mesh_3d_5.cc,v 1.2 2003/10/16 14:18:17 wolf Exp
//    Version: 
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_3d_5.cc  ---------------------------


// check that face orientation flags are properly inherited by
// looking at them and checking their children

#include "../tests.h"
#include "mesh_3d.h"

#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_reordering.h>
#include <grid/grid_generator.h>

#include <fstream>


void check_this (Triangulation<3> &tria)
{
                                   // look at all faces, not only
                                   // active ones
  for (Triangulation<3>::cell_iterator cell=tria.begin();
       cell != tria.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f)
      if (cell->has_children())
        for (unsigned int c=0; c<GeometryInfo<3>::subfaces_per_face; ++c)
          {
            Assert (cell->face_orientation(f) ==
                    cell->child(GeometryInfo<3>::child_cell_on_face(f,c))
                    ->face_orientation(f),
                    ExcInternalError());
            deallog << "Cell << " << cell
                    << ", face " << f
                    << " subface " << c << " is ok."
                    << std::endl;
          }
}


void check (Triangulation<3> &tria)
{
  deallog << "Initial check" << std::endl;
  check_this (tria);
  
  for (unsigned int r=0; r<3; ++r)
    {
      tria.refine_global (1);
      deallog << "Check " << r << std::endl;
      check_this (tria);
    }

  coarsen_global (tria);
  deallog << "Check " << 1 << std::endl;
  check_this (tria);
  
  tria.refine_global (1);
  deallog << "Check " << 2 << std::endl;
  check_this (tria);
}


int main () 
{
  std::ofstream logfile("mesh_3d_5.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  {  
    Triangulation<3> coarse_grid;
    create_two_cubes (coarse_grid);
    check (coarse_grid);
  }
  
  {  
    Triangulation<3> coarse_grid;
    create_L_shape (coarse_grid);
    check (coarse_grid);
  }
  
  {  
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_ball (coarse_grid);
    check (coarse_grid);
  }
  
}

  
  
