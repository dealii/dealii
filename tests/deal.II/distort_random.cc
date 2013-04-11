//----------------------------  distort_random.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  distort_random.cc  ---------------------------

// check GridTools::distort_random


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iomanip>

std::ofstream logfile("distort_random/output");



template <int dim>
void test1 (const bool keep_boundary)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  GridTools::distort_random (0.1, tria, keep_boundary);

  deallog << "dim=" << dim << ", keep_boundary=" << keep_boundary << std::endl;
  GridOut().write_gnuplot (tria, logfile);
}



int main ()
{
  deallog << std::setprecision(4);
  logfile << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test1<1> (true);
  test1<1> (false);
  test1<2> (true);
  test1<2> (false);
  test1<3> (true);
  test1<3> (false);

  return 0;
}

