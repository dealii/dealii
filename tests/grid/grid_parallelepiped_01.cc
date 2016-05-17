// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



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
std::ofstream logfile ("output");

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
  corners[0] = Point<2> (0.5, 0.0);
  corners[1] = Point<2> (0.0, 0.5);

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
