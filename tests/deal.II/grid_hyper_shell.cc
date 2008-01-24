//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


#include "../tests.h"
#include <grid/tria_boundary_lib.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <base/logstream.h>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <iomanip>

std::ofstream logfile("grid_hyper_shell/output");


template<int dim>
void check (double r1, double r2, unsigned int n, bool log)
{
  Point<dim> center;
  Triangulation<dim> tria;
  GridGenerator::hyper_shell (tria, center, r1, r2, n);
  static const HyperShellBoundary<dim> boundary(center);
  tria.set_boundary(0, boundary);

  tria.refine_global(2);
  
  GridOut grid_out;
  if (log)
    grid_out.write_eps (tria, logfile);
  else
    grid_out.write_dx (tria, std::cout);
}


int main()
{
  check<2> (4., 5., 10, true);
  check<3> (3., 5., 6, true);
}
