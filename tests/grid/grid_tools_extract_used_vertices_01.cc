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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Check extract used_vertices.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
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

  auto m = GridTools::extract_used_vertices(tria);

  for (auto &e : m)
    deallog << "Vertex: " << e.first << ": " << e.second << std::endl;
};


int
main()
{
  initlog();
  test<1, 1>();
  test<1, 2>();
  test<1, 3>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
  return 0;
}
