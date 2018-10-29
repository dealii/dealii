// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2017 by the deal.II authors
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
