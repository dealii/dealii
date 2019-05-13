// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

// Check computing bounding boxes of entire triangulations and individual cells

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test_tria_bounding_box()
{
  dealii::Point<dim> p1, p2;
  for (unsigned int k = 0; k < dim; ++k)
    {
      p1[k] = k;
      p2[k] = 2 * k + 1;
    }

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_rectangle(tria, p1, p2);
  tria.refine_global(1);
  const BoundingBox<spacedim> bounding_box =
    GridTools::compute_bounding_box(tria);
  const std::pair<Point<spacedim>, Point<spacedim>> &boundary_points =
    bounding_box.get_boundary_points();

  deallog << boundary_points.first << ", " << boundary_points.second
          << std::endl;
}

int
main(void)
{
  initlog();

  test_tria_bounding_box<1, 1>();
  test_tria_bounding_box<1, 2>();
  test_tria_bounding_box<2, 2>();
  test_tria_bounding_box<1, 3>();
  test_tria_bounding_box<2, 3>();
  test_tria_bounding_box<3, 3>();

  return 0;
}
