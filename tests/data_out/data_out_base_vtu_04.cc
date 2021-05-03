// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

// Test high-order Lagrange VTU output on a 2D shell

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q_generic.h>

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
  // choose some arbitrary radii to reduce possibility of roundoff effects in
  // the double -> float transition and with binary (zlib compressed VTU)
  // output
  GridGenerator::hyper_shell(triangulation,
                             Point<dim>(3.2323343428452032, 2.12432324033),
                             0.53324387343224532,
                             1.032354728342342875235,
                             6);
  triangulation.refine_global(1);

  FE_Q<dim>       fe(cell_order);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double>       vec(dof_handler.n_dofs());
  MappingQGeneric<dim> mapping(cell_order);

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
  data_out.write_vtu(log);
}

int
main()
{
  initlog();
  deallog.get_file_stream() << std::setprecision(9);

  unsigned cell_order = 3;
  check<2>(deallog.get_file_stream(), cell_order);
}
