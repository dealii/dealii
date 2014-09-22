// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



/* Author: Ralf Hartmann, 2005; the MGDoFHandler class had troubles
   when the triangulation had unused vertices before, make sure that
   this is now fixed. */

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/multigrid/mg_dof_handler.h>

#include <fstream>


int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const unsigned int dim=2;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  FE_DGQ<dim> fe(1);
  MGDoFHandler<dim> dof_handler(tria);

  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  Triangulation<dim>::active_cell_iterator
  cell=tria.begin_active(),
  endc=tria.end();
  for (; cell!=endc; ++cell)
    cell->set_coarsen_flag();
  tria.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(fe);

  deallog << "OK" << std::endl;
}
