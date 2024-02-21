// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test high-order Lagrange VTK output on a 2D shell

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <string>

#include "../tests.h"

template <int dim>
void
check(std::ostream &log, unsigned cell_order)
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_shell(triangulation, Point<dim>(), 0.5, 1, 6);
  triangulation.refine_global(1);

  FE_Q<dim>       fe(cell_order);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> vec(dof_handler.n_dofs());
  MappingQ<dim>  mapping(cell_order);

  VectorTools::interpolate(mapping,
                           dof_handler,
                           Functions::SquareFunction<dim>(),
                           vec);

  DataOutBase::VtkFlags flags;
  flags.write_higher_order_cells = true;

  DataOut<dim> data_out;
  data_out.set_flags(flags);
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(vec, "square_function");
  data_out.build_patches(mapping, cell_order, DataOut<dim>::curved_inner_cells);
  data_out.write_vtk(log);
}

int
main()
{
  initlog();

  unsigned cell_order = 3;
  check<2>(deallog.get_file_stream(), cell_order);
}
