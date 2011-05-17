//-----------------------  find_closest_vertex_1.cc  ----------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors and Ralf B. Schulz
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------  find_closest_vertex_1.cc  ----------------------


// take a 3d mesh, take all vertices, shift them a little bit and check that
// we correctly identify the closest vertex position
// The result should be an increasing sequence of numbers

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <fstream>




void check (Triangulation<3> &tria)
{
   const std::vector<Point<3> > &v = tria.get_vertices();
   for(unsigned i=0; i<v.size(); i++)
      deallog << "[" << GridTools::find_closest_vertex(tria, v[i] + Point<3>(0.01, -0.01, 0.01)) << "] ";
         
  deallog << std::endl;
}


int main () 
{
  std::ofstream logfile("find_closest_vertex_1/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Triangulation<3> coarse_grid;
      GridGenerator::hyper_cube (coarse_grid);
      coarse_grid.refine_global (3);
      check (coarse_grid);
    }
  catch (const std::exception &exc)
    {
				       // we shouldn't get here...
      deallog << "Caught an error..." << std::endl;
      deallog << exc.what() << std::endl;
    }
}

  
  
