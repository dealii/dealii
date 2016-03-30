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

std::ofstream logfile ("output");

void check_remove_anisotropy ()
{
  Point<3> (corners) [3];

  corners[0] = Point<3> (1, 0, 0);
  corners[1] = Point<3> (0, 4, 0);
  corners[2] = Point<3> (0, 0, 2);

  const unsigned int n_subdivisions = 1;

  Triangulation<3> triangulation;
  GridGenerator::subdivided_parallelepiped (triangulation, n_subdivisions, corners);
  dealii::GridTools::remove_anisotropy<3>(triangulation,
                                          /*max ratio =*/ 1.2);

  GridOut grid_out;
  grid_out.write_vtk (triangulation, logfile);

  triangulation.clear ();
}

int main ()
{
  check_remove_anisotropy();
}
