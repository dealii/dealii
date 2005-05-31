//----------------------------  rt_1.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors and Oliver Kayser-Herold
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
void test (const unsigned int degree,
           const unsigned int q_order)
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation, -1, 1);

  FE_RaviartThomas<dim> fe (degree);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);

  QGauss<dim-1> q(q_order);
  FEFaceValues<dim> fe_values (fe, q,
                               update_values |
                               update_gradients |
                               update_second_derivatives |
                               update_q_points |
                               update_jacobians);
  fe_values.reinit (dof_handler.begin_active(), 0);

  deallog << "OK" << std::endl;
}


int main () 
{
  std::ofstream logfile("rt_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  for (unsigned int degree=0; degree<3; ++degree)
    for (unsigned int q_order=1; q_order<=3; ++q_order)
      test<2> (degree, q_order);

  return 0;
}
