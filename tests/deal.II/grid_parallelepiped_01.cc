//----------------------------------------------------------------------
//    $Id: grid_parallelepiped.cc 2013-01-12 young $
//    Version: $Name$ 
//
//    Copyright (C) 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


#include "../tests.h"
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <iomanip>

// Output
std::ofstream logfile ("grid_parallelepiped_01/output");

// As sketched in the deal.II docs, the parallelepiped class is just a
// hyper_rectangle in 1d and a parallelogram in 2d. That can checked
// by simple comparison to a reference triangulation.

// Here is the implementation in 1d:
void check_1d_parallelepiped_by_comparison (bool log)
{
  // Data structure defining dim coordinates that make up a
  // parallelepiped.
  Point<1> (corners) [1];
  corners[0] = Point<1> (0.5);
  
  Triangulation<1> triangulation_parallelepiped;
  GridGenerator::parallelepiped (triangulation_parallelepiped, corners, false);
  
  Triangulation<1> triangulation_cube;
  GridGenerator::hyper_cube (triangulation_cube, 0., 0.5);

  if (log)
    {
      logfile << "\ncheck 1d parallelepiped (hyper_cube): ";
      if (GridTools::have_same_coarse_mesh (triangulation_parallelepiped,
					    triangulation_cube))
	logfile << "OK";
      else
	logfile << "not OK... coarse grids are different but they should be the same";
    }
}


// Here is the implementation in 2d:
void check_2d_parallelepiped_by_comparison (bool log)
{
  // build corners for this particular dim that are known to give the
  // same output order as parallelogram:
  Point<2> (corners) [2];
  corners[0] = Point<2> (0.0, 0.5);
  corners[1] = Point<2> (0.5, 0.0);
  
  Triangulation<2> triangulation_parallelepiped;
  GridGenerator::parallelepiped (triangulation_parallelepiped, corners, false);

  Triangulation<2> triangulation_parallelogram;
  GridGenerator::parallelogram (triangulation_parallelogram, corners, false);

  if (log)
    {
      logfile << "\ncheck 2d parallelepiped (parallelogram): ";
      if (GridTools::have_same_coarse_mesh (triangulation_parallelepiped,
					    triangulation_parallelogram))
	logfile << "OK";

      else
	logfile << "not OK... coarse grids are different but they should be the same";
    }
}

int main ()
{
  // Check parallelepiped 
  check_1d_parallelepiped_by_comparison (true);
  check_2d_parallelepiped_by_comparison (true);
  logfile << "\n";
}
