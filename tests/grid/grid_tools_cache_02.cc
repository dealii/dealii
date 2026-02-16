// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Validate grid_tools_cache. Different construction order. Check signal.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  deallog << "dim = " << dim << ", spacedim = " << spacedim << std::endl;

  Triangulation<dim>    tria;
  GridTools::Cache<dim> cache(tria);

  GridGenerator::hyper_cube(tria);

  {
    auto m = cache.get_vertex_to_cell_map();

    unsigned int i = 0;
    for (auto &v : m)
      {
        for (auto &c : v)
          deallog << "Vertex " << i << ", cell " << c << std::endl;
        ++i;
      }
  }
  deallog << "After one refinement. " << std::endl;

  tria.refine_global(1);
  {
    auto m = cache.get_vertex_to_cell_map();

    unsigned int i = 0;
    for (auto &v : m)
      {
        for (auto &c : v)
          deallog << "Vertex " << i << ", cell " << c << std::endl;
        ++i;
      }
  }
}


int
main()
{
  initlog();

  test<2, 2>();

  return 0;
}
