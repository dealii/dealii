// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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


// when we coarsen away a whole level, the level data structures
// remain but contains only unused cells. make sure that tria.last()
// and tria.last_active() still produce something sensible in that
// case

#include "../tests.h"
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/logstream.h>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iomanip>
#include <cstdio>

std::ofstream logfile("output");



template <int dim>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);
  for (unsigned int i=0; i<2; ++i)
    {
      for (typename Triangulation<dim>::active_cell_iterator
	     cell = tria.begin_active(2); cell != tria.end(); ++cell)
	cell->set_coarsen_flag ();
      tria.execute_coarsening_and_refinement ();
    }

  deallog << tria.n_levels() << ' ' << tria.n_global_levels() << ' '
	  << tria.last() << ' '
	  << tria.last_active()
	  << std::endl;
}


int main ()
{
  deallog << std::setprecision(2);
  logfile << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1>();
  test<2>();
  test<3>();

  return 0;
}
