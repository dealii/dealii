//----------------------------  crash_02.cc  ---------------------------
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
//----------------------------  crash_02.cc  ---------------------------


// check a crash in hp::DoFHandler because the
// DoFHandler::active_fe_indices array isn't initialized when
// attaching to a triangulation


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_dgq.h>

#include <fstream>


template <int dim>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back (FE_DGQ<dim> (1));

  hp::DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs (fe_collection);
  
  dof_handler.begin_active()->set_active_fe_index(0);
}



int main ()
{
  std::ofstream logfile("crash_02/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  test<1> ();
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
