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

void check_remove_hanging_nodes ()
{
  Point<2> corners[2];

  corners[0] = Point<2> (1, 0);
  corners[1] = Point<2> (0, 4);
  const unsigned int n_subdivisions = 1;


  Triangulation<2> tria;
  GridGenerator::subdivided_parallelepiped (tria, n_subdivisions, corners);

  tria.refine_global();

  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  tria.refine_global();

  dealii::GridTools::remove_hanging_nodes(tria, /*isotropic=*/false);

  GridOut grid_out;
  grid_out.write_vtk (tria, logfile);

  tria.clear ();
}

int main ()
{
  check_remove_hanging_nodes();
}
