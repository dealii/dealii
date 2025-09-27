// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Output a field constant in time

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out_stack.h>

#include "../tests.h"



template <int dim>
void
run()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(1);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  // create a continuous field over
  // this DoFHandler
  Vector<double> v(dof_handler.n_dofs());
  v               = 1.;
  v(v.size() / 2) = 2.;

  // output this field using
  // DataOutStack. the result should
  // be a continuous field again
  DataOutStack<dim> data_out_stack;
  data_out_stack.declare_data_vector("solution", DataOutStack<dim>::dof_vector);
  data_out_stack.new_parameter_value(1, 1);
  data_out_stack.attach_dof_handler(dof_handler);
  data_out_stack.add_data_vector(v, "solution");
  data_out_stack.build_patches(2);
  data_out_stack.finish_parameter_value();

  data_out_stack.write_dx(deallog.get_file_stream());
}



int
main()
{
  initlog();

  run<1>();
  run<2>();
  return 0;
}
