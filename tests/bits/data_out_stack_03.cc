// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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


// Output a field constant in time

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out_stack.h>
#include <fstream>
#include <iomanip>

#include <deal.II/base/logstream.h>




template <int dim>
void run ()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (1);

  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs (fe);

  // create a continuous field over
  // this DoFHandler
  Vector<double> v(dof_handler.n_dofs());
  v = 1.;
  v(v.size()/2) = 2.;

  // output this field using
  // DataOutStack. the result should
  // be a continuous field again
  DataOutStack<dim> data_out_stack;
  data_out_stack.declare_data_vector ("solution",
                                      DataOutStack<dim>::dof_vector);
  data_out_stack.new_parameter_value (1,1);
  data_out_stack.attach_dof_handler (dof_handler);
  data_out_stack.add_data_vector (v, "solution");
  data_out_stack.build_patches (2);
  data_out_stack.finish_parameter_value ();

  data_out_stack.write_dx (deallog.get_file_stream());
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  run<1> ();
  run<2> ();
  return 0;
}
