// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// this test failed with an internal error somewhere in the coarsening
// functions

#include "../tests.h"
#include "../grid/mesh_3d.h"

#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_generator.h>

#include <fstream>


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Triangulation<3> coarse_grid;
  create_L_shape (coarse_grid);

  // refine once, then unrefine again
  coarse_grid.refine_global (1);
  for (Triangulation<3>::active_cell_iterator c=coarse_grid.begin_active();
       c != coarse_grid.end(); ++c)
    c->set_coarsen_flag ();
  coarse_grid.execute_coarsening_and_refinement ();

  deallog << "ok." << std::endl;
}



