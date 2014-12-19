// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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

// The simplest test case is to create a parallelepiped grid with a
// number of subdivisions and output the result.
template<int dim>
void check_subdivided_parallelepiped (bool colorize, bool log)
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
      corners[0] = Point<dim> (0.25, 0.50);
      corners[1] = Point<dim> (0.50, 0.25);
      break;

    case 3:
      corners[0] = Point<dim> (0.25, 0.50, 0.50);
      corners[1] = Point<dim> (0.50, 0.25, 0.50);
      corners[2] = Point<dim> (0.50, 0.50, 0.25);
      break;

    default:
      Assert (false, ExcInternalError ());
    }

  // The number of subdivisions can be anything reasonable:
  const unsigned int n_subdivisions = (2*dim+1);

  Triangulation<dim> triangulation;
  GridGenerator::subdivided_parallelepiped (triangulation, n_subdivisions, corners, colorize);

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
  check_subdivided_parallelepiped<1> (false, true);
  check_subdivided_parallelepiped<2> (false, true);
  check_subdivided_parallelepiped<3> (true,  true);
}
