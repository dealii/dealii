//----------------------------  rt_1.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors and
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  rt_1.cc  ---------------------------


// there was a bug in the RT element that Oliver Kayser-Herold fixed in
// January 2005. Check this

#include "../tests.h"
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_raviart_thomas.h>	
#include <fe/fe_values.h>
#include <iostream>
#include <fstream>


template <int dim>
void test ()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation, -1, 1);

  DoFHandler<dim> dof_handler (triangulation);
  FE_RaviartThomas<dim> fe (1);
  dof_handler.distribute_dofs (fe);

  QGauss<dim-1> q(2);
  FEFaceValues<dim> fe_values (fe, q, update_values|update_gradients);
  fe_values.reinit (dof_handler.begin_active(), 0);

  deallog << "OK" << std::endl;
}


int main () 
{
  std::ofstream logfile("rt_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<2> ();
  test<3> ();
  return 0;
}
