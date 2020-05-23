// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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



// Test for a patch which makes sure that merge_triangulations() does
// not forget about where the boundary is (the problem has
// actually been in create_triangulation(), which assigned
// numbers::internal_face_boundary_id to boundary faces)

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"


int
main()
{
  initlog();

  const unsigned int spacedim = 3;

  Triangulation<spacedim> tria_0, tria_1, tria_2;

  GridGenerator::hyper_cube(tria_0);
  GridGenerator::hyper_cube(tria_1);
  Tensor<1, spacedim> shift;
  shift[0] = 1.0;
  GridTools::shift(shift, tria_1);

  tria_0.set_all_manifold_ids(0);
  tria_1.set_all_manifold_ids(0);
  GridGenerator::merge_triangulations(tria_0, tria_1, tria_2, 1e-12, true);

  unsigned int boundary_face_count = 0;
  for (const auto &cell : tria_2.active_cell_iterators())
    for (const unsigned int f : GeometryInfo<spacedim>::face_indices())
      if (cell->face(f)->at_boundary())
        ++boundary_face_count;

  if (boundary_face_count != 10)
    deallog << "Found " << boundary_face_count
            << " boundary faces. However, there should be 10 of them!"
            << std::endl;
  else
    deallog << "OK!" << std::endl;
}
