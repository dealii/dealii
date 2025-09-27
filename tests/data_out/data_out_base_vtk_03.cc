// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
  check<1>(deallog.get_file_stream(), cell_order);
  check<2>(deallog.get_file_stream(), cell_order);
  check<3>(deallog.get_file_stream(), cell_order);
}
