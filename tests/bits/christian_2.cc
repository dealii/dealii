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


// check one aspect of DataOutStack in 2+1d

#include "../tests.h"
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_out_stack.h>
#include <fstream>
#include <iomanip>


int main()
{
  std::ofstream logfile("output");
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
