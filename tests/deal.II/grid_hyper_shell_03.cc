//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

// generate and refine a hyper shell in 3d with 12 cells. The other tests of
// this series test the other possibilities for the number of cells

#include "../tests.h"
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/logstream.h>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <iomanip>

std::ofstream logfile("grid_hyper_shell_03/output");


template<int dim>
void check (double r1, double r2, unsigned int n)
{
  Point<dim> center;
  Triangulation<dim> tria (Triangulation<dim>::none);
  GridGenerator::hyper_shell (tria, center, r1, r2, n);
  static const HyperShellBoundary<dim> boundary(center);
  tria.set_boundary(0, boundary);
  tria.refine_global(1);

  GridOut grid_out;
  grid_out.write_gnuplot (tria, logfile);
}


int main()
{
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  check<3> (.5, 1, 12);
}
