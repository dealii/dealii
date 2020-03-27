// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


// Test GridTools::partition_triangulation_zorder

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test(const int n_refinements, const int n_partitions, const bool blocked)
{
  // create serial triangulation
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_refinements);

  GridTools::partition_triangulation_zorder(n_partitions, tria, blocked);

  for (const auto &cell : tria.active_cell_iterators())
    deallog << cell->subdomain_id() << " ";
  deallog << std::endl;
}

int
main(int argc, char *argv[])
{
  initlog();

  if (true)
    {
      deallog.push("2d");
      test<2>(1, 4, false);
      deallog.pop();
    }
  if (true)
    {
      deallog.push("2d-pdt");
      test<2>(1, 4, true);
      deallog.pop();
    }
  if (true)
    {
      deallog.push("3d");
      test<3>(1, 8, false);
      deallog.pop();
    }
  if (true)
    {
      deallog.push("3d-pdt");
      test<3>(1, 8, true);
      deallog.pop();
    }
}
