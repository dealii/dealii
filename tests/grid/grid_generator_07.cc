// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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

// Test GridGenerator::create_triangulation_with_removed_cells

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

template <int dim>
void
test(std::ostream& out)
{
  Triangulation<dim> triangulation;
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  // remove all cells but the first. this is the hardest case to handle as it
  // makes a bunch of vertices unused
  std::set<typename Triangulation<dim>::active_cell_iterator> cells_to_remove;
  for(typename Triangulation<dim>::active_cell_iterator cell
      = ++triangulation.begin_active();
      cell != triangulation.end();
      ++cell)
    cells_to_remove.insert(cell);

  GridGenerator::create_triangulation_with_removed_cells(
    triangulation, cells_to_remove, tr);
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
