//----------------------------  mesh_3d_4.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_3d_4.cc  ---------------------------


// check that face orientation flags are properly inherited by
// counting them

#include "mesh_3d.h"

#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_reordering.h>
#include <grid/grid_generator.h>

#include <fstream>


unsigned int count_wrong_faces (const Triangulation<3> &tria)
{
  unsigned int count = 0;
  
                                   // count faces with "wrong"
                                   // orientation
  for (Triangulation<3>::active_cell_iterator cell=tria.begin_active();
       cell != tria.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f)
      if (cell->face_orientation(f) == false)
        ++count;
  return count;
}

  

void check (Triangulation<3> &tria)
{
  const unsigned int initial_count = count_wrong_faces (tria);
  for (unsigned int r=0; r<3; ++r)
    {
      tria.refine_global (1);
      const unsigned int count = count_wrong_faces (tria);
      deallog << "'Wrong' faces = " << count << std::endl;
      Assert (count == initial_count * (4<<(2*r)),
              ExcInternalError());
    }
  
  {
    coarsen_global (tria);
    const unsigned int count = count_wrong_faces (tria);
    deallog << "'Wrong' faces = " << count << std::endl;
    Assert (count == initial_count * (4<<(2*1)),
            ExcInternalError());
  }
  
  {
    tria.refine_global (1);
    const unsigned int count = count_wrong_faces (tria);
    deallog << "'Wrong' faces = " << count << std::endl;
    Assert (count == initial_count * (4<<(2*2)),
            ExcInternalError());
  }
}

  

int main () 
{
  std::ofstream logfile("mesh_3d_4.output");
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

  
  
