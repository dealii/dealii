// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test GridGenerator::create_triangulation_with_removed_cells

#include <deal.II/base/tensor.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test(std::ostream &out)
{
  Triangulation<dim> triangulation;
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  // remove all cells but the first. this is the hardest case to handle as it
  // makes a bunch of vertices unused
  std::set<typename Triangulation<dim>::active_cell_iterator> cells_to_remove;
  for (typename Triangulation<dim>::active_cell_iterator cell =
         std::next(triangulation.begin_active());
       cell != triangulation.end();
       ++cell)
    cells_to_remove.insert(cell);

  GridGenerator::create_triangulation_with_removed_cells(triangulation,
                                                         cells_to_remove,
                                                         tr);
  GridOut go;
  go.write_gnuplot(tr, out);
}


int
main()
{
  initlog();

  test<2>(deallog.get_file_stream());
  test<3>(deallog.get_file_stream());
}
