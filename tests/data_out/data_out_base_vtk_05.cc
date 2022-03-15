// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Test high-order Lagrange VTK output on a 3D shell

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
  GridGenerator::hyper_shell(triangulation, Point<dim>(1, 0, 0), 0.5, 1, 6);

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
  check<3>(deallog.get_file_stream(), cell_order);
}
