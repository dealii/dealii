// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
    deallog << cell->subdomain_id() << ' ';
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
