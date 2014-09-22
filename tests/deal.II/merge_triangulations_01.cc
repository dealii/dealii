// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/logstream.h>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iomanip>
#include <cstdio>

std::ofstream logfile("output");



// testcase=0:
// create two cubes; translate them so that no vertices overlap
//
// testcase=1:
// create two cubes; translate them so that a whole face overlaps
//
// testcase=2:
// create two cubes; translate them so that exactly one vertices overlaps
template <int dim>
void test (const int testcase)
{
  Triangulation<dim> tria_1, tria_2, tria_3;
  GridGenerator::hyper_cube(tria_1);
  GridGenerator::hyper_cube(tria_2);
  Point<dim> shift;
  switch (testcase)
    {
    case 0:
      shift[0] = 2;
      break;
    case 1:
      shift[0] = 1;
      break;
    case 2:
      for (unsigned int d=0; d<dim; ++d)
        shift[d] = 1;
      break;
    default:
      Assert (false, ExcNotImplemented());
    }
  GridTools::shift (shift, tria_2);

  // fill tria_3 with something, to
  // make sure that the function we
  // call later can deal with prior
  // content
  GridGenerator::hyper_cube(tria_3);

  // now merge triangulations
  GridGenerator::merge_triangulations (tria_1, tria_2, tria_3);

  GridOut().write_gnuplot (tria_3, logfile);

  deallog << "     Total number of cells        = " << tria_3.n_cells() << std::endl
          << "     Total number of vertices = " << tria_3.n_used_vertices() << std::endl;
}


int main ()
{
  deallog << std::setprecision(2);
  logfile << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  for (unsigned int t=0; t<3; ++t)
    {
      test<2> (t);
      test<3> (t);
    }

  return 0;
}
