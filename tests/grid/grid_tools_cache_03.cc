// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check extract used_vertices and find_closest_vertex using a cache

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test(const Point<spacedim> &p)
{
  deallog << "dim: " << dim << ", spacedim: " << spacedim << std::endl;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);

  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  for (auto cell = tria.begin_active(2); cell != tria.end(2); ++cell)
    cell->set_coarsen_flag();

  tria.execute_coarsening_and_refinement();

  GridTools::Cache<dim, spacedim> cache(tria);

  auto m = cache.get_used_vertices();

  for (auto &e : m)
    deallog << "Vertex: " << e.first << ": " << e.second << std::endl;

  auto i = GridTools::find_closest_vertex(m, p);
  deallog << "Closest vertex to " << p << ", v[" << i << "] :" << m[i]
          << std::endl;
};


int
main()
{
  initlog();
  test<2, 2>(Point<2>(.2, .2));
  test<2, 2>(Point<2>(.6, .9));
  return 0;
}
