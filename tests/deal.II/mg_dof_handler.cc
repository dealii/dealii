//----------------------------  mg_dof_handler.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mg_dof_handler.cc  ---------------------------


/* Author: Ralf Hartmann, 2005; the MGDoFHandler class had troubles
   when the triangulation had unused vertices before, make sure that
   this is now fixed. */

#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_out.h>
#include <fe/fe_dgq.h>
#include <multigrid/mg_dof_handler.h>

#include <fstream>


int main ()
{
  std::ofstream logfile("mg_dof_handler.output");
  logfile.precision(2);
  
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
