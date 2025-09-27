// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check one aspect of DataOutStack in 2+1d

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out_stack.h>

#include "../tests.h"


int
main()
{
  initlog();

  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria, 0, 1);
  FE_Q<2> fe(1);

  DoFHandler<2> dof(tria);
  dof.distribute_dofs(fe);

  Vector<double> sol(dof.n_dofs());
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    sol(i) = i;

  // test output using DataOutStack
  DataOutStack<2> data_out_stack;
  data_out_stack.declare_data_vector("dof", DataOutStack<2>::dof_vector);
  data_out_stack.new_parameter_value(2.5, 1.);
  data_out_stack.attach_dof_handler(dof);
  data_out_stack.add_data_vector(sol, "dof");
  data_out_stack.build_patches();
  data_out_stack.finish_parameter_value();

  data_out_stack.write_gnuplot(deallog.get_file_stream());
}
