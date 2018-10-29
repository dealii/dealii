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

// Test high-order Lagrange VTK output on a unit cube.
// Added by Alexander Grayver

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>

#include <string>

#include "../tests.h"

template <int dim>
void
check(std::ostream &log, unsigned cell_order)
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);

  DataOutBase::VtkFlags flags;
  flags.write_higher_order_cells = true;

  DataOut<dim> data_out;
  data_out.set_flags(flags);
  data_out.attach_triangulation(triangulation);
  data_out.build_patches(cell_order);
  data_out.write_vtk(log);
}

int
main()
{
  initlog();

  unsigned cell_order = 4;
  check<2>(deallog.get_file_stream(), cell_order);
  check<3>(deallog.get_file_stream(), cell_order);
}
