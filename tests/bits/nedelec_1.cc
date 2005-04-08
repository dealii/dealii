//----------------------------  nedelec_1.cc  ---------------------------
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
//----------------------------  nedelec_1.cc  ---------------------------


// this test is modeled after after the nedelec_1 test, since that one failed. I
// just wanted to check that the nedelec element that has much of the same
// code also works. turns out, all was fine

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
#include <fe/fe_nedelec.h>	
#include <fe/fe_values.h>
#include <iostream>
#include <fstream>


template <int dim>
void test ()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation, -1, 1);

  FE_Nedelec<dim> fe (1);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);

  QGauss<dim-1> q(2);
  FEFaceValues<dim> fe_values (fe, q, update_values|update_gradients);
  fe_values.reinit (dof_handler.begin_active(), 0);

  deallog << "OK" << std::endl;
}


int main () 
{
  std::ofstream logfile("nedelec_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2> ();
  test<3> ();
  
  return 0;
}
