// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Generate a grid, refine it once, flatten it and output the result.
//
// This test checks that we are also copying the boundary and manifold
// ids of edges in 3d.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test()
{
  deallog << dim << 'd' << std::endl;

  // Create a triangulation, then set the boundary and manifold ids of
  // boundary faces and, in 3d, of boundary edges
  Triangulation<dim> tria1;
  GridGenerator::hyper_cube(tria1);
  tria1.refine_global(1);

  for (const auto &cell : tria1.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary())
        {
          face->set_boundary_id(101);
          face->set_manifold_id(102);

          if (dim == 3)
            for (unsigned int l = 0; l < face->n_lines(); ++l)
              {
                face->line(l)->set_boundary_id(201);
                face->line(l)->set_manifold_id(202);
              }
        }


  // Then flatten the triangulation, and output the boundary and
  // manifold ids of all objects on the boundary
  Triangulation<dim> tria2;
  GridGenerator::flatten_triangulation(tria1, tria2);

  for (const auto &cell : tria2.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary())
        {
          deallog << "Face " << face << std::endl;
          deallog << "  center = " << face->center() << std::endl;
          deallog << "  boundary id = " << face->boundary_id() << std::endl;
          deallog << "  manifold id = " << face->manifold_id() << std::endl;

          if (dim == 3)
            for (unsigned int l = 0; l < face->n_lines(); ++l)
              {
                deallog << "  Edge " << face->line(l) << std::endl;
                deallog << "    boundary id = " << face->line(l)->boundary_id()
                        << std::endl;
                deallog << "    manifold id = " << face->line(l)->manifold_id()
                        << std::endl;
              }
        }

  deallog << std::endl << std::endl;
}

int
main()
{
  initlog();

  test<2>();
  test<3>();
}
