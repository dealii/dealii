// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

// Test GridGenerator::concentric_hyper_shell with skewness = 2.0

int
main()
{
  initlog();

  // 2D test
  {
    Triangulation<2> triangulation;
    const Point<2>   center(0.0, 1.0);
    GridGenerator::concentric_hyper_shells(
      triangulation, center, 1.0, 4.0, 3u, 2.0, 4, true);

    for (const auto &cell : triangulation.active_cell_iterators())
      {
        deallog << "vertices: " << cell->vertex(0) << ", " << cell->vertex(1)
                << ", " << cell->vertex(2) << ", " << cell->vertex(3)
                << std::endl;

        bool manifold_ids_are_zero = cell->manifold_id() == 0;
        for (const unsigned int face_n : GeometryInfo<2>::face_indices())
          manifold_ids_are_zero &= cell->face(face_n)->manifold_id() == 0;
        AssertThrow(manifold_ids_are_zero, ExcInternalError());

        for (const unsigned int face_n : GeometryInfo<2>::face_indices())
          if (cell->face(face_n)->at_boundary())
            deallog << "boundary face center distance to origin: "
                    << (cell->face(face_n)->center(/*respect_manifold*/ true) -
                        center)
                         .norm()
                    << " boundary id: " << cell->face(face_n)->boundary_id()
                    << std::endl;
      }
  }

  // 3D test
  {
    Triangulation<3> triangulation;
    const Point<3>   center(0.0, 1.0, 2.0);
    GridGenerator::concentric_hyper_shells(
      triangulation, center, 1.0, 5.0, 2u, 2.0, 0, true);

    for (const auto &cell : triangulation.active_cell_iterators())
      {
        deallog << "vertices: " << cell->vertex(0) << ", " << cell->vertex(1)
                << ", " << cell->vertex(2) << ", " << cell->vertex(3) << ", "
                << cell->vertex(4) << ", " << cell->vertex(5) << ", "
                << cell->vertex(6) << ", " << cell->vertex(7) << std::endl;

        bool manifold_ids_are_zero = cell->manifold_id() == 0;
        for (const unsigned int face_n : GeometryInfo<3>::face_indices())
          manifold_ids_are_zero &= cell->face(face_n)->manifold_id() == 0;
        AssertThrow(manifold_ids_are_zero, ExcInternalError());

        for (const unsigned int face_n : GeometryInfo<3>::face_indices())
          if (cell->face(face_n)->at_boundary())
            deallog << "boundary face center distance to origin: "
                    << (cell->face(face_n)->center(/*respect_manifold*/ true) -
                        center)
                         .norm()
                    << " boundary id: " << cell->face(face_n)->boundary_id()
                    << std::endl;
      }
  }
}
