// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

// Check if Triangulation<1>::get_boundary_indicators() works for 1d grids.


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  Triangulation<1>   triangulation;
  GridGenerator::hyper_cube (triangulation, -1, 1);
  const std::vector<types::boundary_id> indicators = triangulation.get_boundary_indicators();
  for (unsigned int i=0; i<indicators.size(); ++i)
    deallog << int (indicators[i]) << std::endl;

  return 0;
}
