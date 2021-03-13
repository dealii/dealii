// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

// Test assign_boundary_ids by setting ids to two faces.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  // assign b-id of 1 to (x == 0) and b-id of 2 to (x == 1)
  // Two lambdas are passed. One of them accepts a Point<spacedim>
  // however, the other one can accept anything. The function should
  // work with either of these lambda definitions.
  GridTools::assign_boundary_ids(
    tria,
    {{1, [](const Point<spacedim> &p) { return p[0] == 0.; }},
     {2, [](const auto &p) { return p[0] == 1.; }}});
  deallog << "hyper_cube<" << dim << ", " << spacedim << ">" << std::endl;
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      for (auto f : cell->face_indices())
        if (cell->face(f)->at_boundary())
          deallog << "face " << f << " B-id " << cell->face(f)->boundary_id()
                  << std::endl;
}

int
main()
{
  initlog();

  test<1, 1>();
  test<1, 2>();
  test<2, 2>();
  test<1, 3>();
  test<2, 3>();
  test<3, 3>();
  deallog << "OK" << std::endl;
  return 0;
}
