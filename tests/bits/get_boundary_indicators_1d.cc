//----------------------------  get_boundary_indicators_1d.cc  ---------------------------
//    $Id: testsuite.html 23951 2011-07-20 11:44:09Z bangerth $
//
//    Copyright (C) 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  get_boundary_indicators_1d.cc  ---------------------------


#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

// Check if Triangulation<1>::get_boundary_indicators() works for 1d grids.


int main ()
{
  std::ofstream logfile("get_boundary_indicators_1d/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  Triangulation<1>   triangulation;
  GridGenerator::hyper_cube (triangulation, -1, 1);
  const std::vector<unsigned char> indicators = triangulation.get_boundary_indicators();
  for (unsigned int i=0; i<indicators.size(); ++i)
    deallog << int (indicators[i]) << std::endl;

  return 0;
}
