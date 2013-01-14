//----------------------------------------------------------------------
//    $Id: grid_parallelepiped.cc 2013-01-07 young $
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
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <iomanip>

// Output
std::ofstream logfile ("grid_parallelepiped/output");

// The simplest test case is to create a parallelepiped grid, output
// the result, and hope for the best.
template<int dim>
void check_parallelepiped (bool colorize, bool log)
{
  // Data structure defining dim coordinates that make up a
  // parallelepiped.
  Point<dim> (corners) [dim];

  // build corners for this particular dim:
  switch (dim)
    {
    case 1:
      corners[0] = Point<dim> (0.5);
      break;

    case 2:
      corners[0] = Point<dim> (0.0, 0.5);
      corners[1] = Point<dim> (0.5, 0.0);
      break;
      
    case 3:
      corners[0] = Point<dim> (0.0, 0.5, 0.5);
      corners[1] = Point<dim> (0.5, 0.0, 0.5);
      corners[2] = Point<dim> (0.5, 0.5, 0.0);
      break;

    default:
      Assert (false, ExcInternalError ());
    }
  
  Triangulation<dim> triangulation;
  GridGenerator::parallelepiped (triangulation, corners, colorize);
  
  GridOut grid_out;

  if (log)
    grid_out.write_gnuplot (triangulation, logfile);

  else
    grid_out.write_gnuplot (triangulation, std::cout);

  triangulation.clear ();
}

int main ()
{
  // Check parallelepiped 
  check_parallelepiped<1> (false, true);
  check_parallelepiped<2> (false, true);
  check_parallelepiped<3> (true,  true);
}
