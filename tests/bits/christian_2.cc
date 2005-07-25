//----------------------------  christian_2.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors and Christian Kamm
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  christian_2.cc  ---------------------------

// check one aspect of DataOutStack in 2+1d

#include "../tests.h"
#include <lac/vector.h>
#include <grid/grid_generator.h>
#include <fe/fe_q.h>
#include <dofs/dof_handler.h>
#include <numerics/data_out_stack.h>
#include <fstream>
#include <iostream>


int main()
{
  std::ofstream logfile("christian_2.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria, 0, 1);
  FE_Q<2> fe(1);
 
  DoFHandler<2> dof(tria);
  dof.distribute_dofs(fe);

  Vector<double> sol(dof.n_dofs());
  for (unsigned int i=0; i<dof.n_dofs(); ++i)
    sol(i) = i;

                                   // test output using DataOutStack
  DataOutStack<2> data_out_stack;
  data_out_stack.declare_data_vector("dof", DataOutStack<2>::dof_vector);
  data_out_stack.new_parameter_value(2.5,1.);
  data_out_stack.attach_dof_handler(dof);
  data_out_stack.add_data_vector(sol, "dof");
  data_out_stack.build_patches();
  data_out_stack.finish_parameter_value();

  data_out_stack.write_gnuplot(logfile);
}
