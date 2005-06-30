//----------------------------  point_inside_1.cc  ---------------------------
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


template <int dim>
void check ()
{
  Triangulation<dim> triangulation;
  
  GridGenerator::hyper_cube (triangulation);
  
				   // Now get the cell
  const Triangulation<dim>::cell_iterator cell = triangulation.begin();
    
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
    
  const bool expected2d[] = {1,0,0,1,0,0,1,1,1,1,0};
  const bool expected3d[] = {1,0,0,0,0,0,0,1,1,1,0};
  const bool *expected=dim==2 ? expected2d : expected3d;    
  for (int i=0; i<11; i++)
    {
      Point<dim> testpoint;
      testpoint(0)=testcoord[i][0];
      testpoint(1)=testcoord[i][1];
      if (dim==3)
	testpoint(2)=testcoord[i][2];

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
  deallog.threshold_double(1.e-10);

  check<2> ();
  check<3> ();
}
