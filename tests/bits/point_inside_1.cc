//----------------------------  point_inside_1.cc  ---------------------------
//    point_inside_1.cc,v 1.2 2003/09/24 15:24:58 wolf Exp
//    Version: 
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  point_inside_1.cc  ---------------------------


// check TriaAccessor<3>::point_inside for a cell that is aligned with
// the coordinate axes
//
// this program is a modified version of one by Joerg Weimar,
// TU Braunschweig

#include "../tests.h"
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <fstream>


void check ()
{
  Triangulation<3> triangulation;
  
  GridGenerator::hyper_cube (triangulation);
  
				   // Now get the cell
  const Triangulation<3>::cell_iterator cell = triangulation.begin();
    
  double testcoord[11][3] = {{0.5,0.5,0.5},
			     {2,0.5,0.5},
			     {0.5,2,0.5},
			     {0.5,0.5,2},
			     {-2,0.5,0.5},
			     {0.5,-2,0.5},
			     {0.5,0.5,-2},
			     {0.9,0.9,0.9},
			     {1.0,0.5,0.5},
			     {0.9999999,0.5,0.5},
			     {1.0000001,0.5,0.5}  };
    
  int expected[] = {1,0,0,0,0,0,0,1,1,1,0};
    
  for (int i=0; i<11; i++)
    {
      const Point<3> testpoint(testcoord[i][0],
			       testcoord[i][1],
			       testcoord[i][2]);
      bool res = cell->point_inside(testpoint);
      deallog << testpoint << " inside " << res <<std::endl;
      Assert (res == expected[i], ExcInternalError());
    }
}


int main () 
{
  std::ofstream logfile("point_inside_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  check ();
}
