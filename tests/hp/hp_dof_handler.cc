//----------------------------  hp_dof_handler.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_dof_handler.cc  ---------------------------


/* Author: Ralf Hartmann, 2005, O. Kayser-Herold, simply modified
  the mg_dof_handler.cc test for the hp::DoFHandler. */

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/hp/dof_handler.h>

#include <fstream>


int main ()
{
  std::ofstream logfile("hp_dof_handler/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const unsigned int dim=2;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back (FE_DGQ<dim> (1));

  hp::DoFHandler<dim> dof_handler(tria);
  
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  Triangulation<dim>::active_cell_iterator
    cell=tria.begin_active(),
    endc=tria.end();
  for (; cell!=endc; ++cell)
    cell->set_coarsen_flag();
  tria.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(fe_collection);

  deallog << "OK" << std::endl;
}
