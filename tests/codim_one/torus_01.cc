//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2008, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// create a torus mesh and refine it.

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>

// all include files you need here

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <fstream>
#include <string>

std::ofstream logfile("torus_01/output");


int main ()
{
  const int dim = 2;
  const int spacedim = 3;

  deallog.attach(logfile);
  deallog.depth_console(0);

  TorusBoundary<dim, spacedim> boundary (1.5, .5);
  Triangulation<dim, spacedim> tria;
  tria.set_boundary (0, boundary);

  GridGenerator::torus (tria, 1.5, .5);
  tria.refine_global(2);

  GridOut grid_out;
  grid_out.write_gnuplot (tria, logfile);

  return 0;
}

