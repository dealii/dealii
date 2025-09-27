// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Validate grid_tools_cache

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  deallog << "dim = " << dim << ", spacedim = " << spacedim << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  GridTools::Cache<dim> cache(tria);

  auto m = cache.get_vertex_to_cell_map();

  unsigned int i = 0;
  for (auto &v : m)
    {
      for (auto &c : v)
        deallog << "Vertex " << i << ", cell " << c << std::endl;
      ++i;
    }
}


int
main()
{
  initlog();

  test<1, 1>();
  test<1, 2>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();

  return 0;
}
