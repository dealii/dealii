/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2021 - 2022 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Like convert_hypercube_to_simplex_mesh_01, but also refines the grid once.
 */

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <fstream>
#include <iostream>

#include "../tests.h"


template <int dim, int spacedim>
void
check()
{
  Triangulation<dim, spacedim> in_tria;
  GridGenerator::quarter_hyper_ball(in_tria);
  in_tria.refine_global(1);

  // Triangulation<dim, spacedim> in_tria_x;
  // GridGenerator::flatten_triangulation(in_tria, in_tria_x);

  Triangulation<dim, spacedim> out_tria;
  GridGenerator::convert_hypercube_to_simplex_mesh(in_tria, out_tria);

  for (const auto &cell : out_tria.active_cell_iterators())
    {
      deallog << "cell=" << cell << ", material_id=" << cell->material_id()
              << ", manifold_id=" << cell->manifold_id() << std::endl;
      for (const auto &face : cell->face_iterators())
        {
          deallog << "  face=" << face
                  << ", boundary_id=" << face->boundary_id()
                  << ", manifold_id=" << face->manifold_id() << std::endl;
          for (unsigned int l = 0; l < face->n_lines(); ++l)
            {
              const auto edge = face->line(l);
              deallog << "    edge=" << edge
                      << ", boundary_id=" << edge->boundary_id()
                      << ", manifold_id=" << edge->manifold_id() << std::endl;
            }
        }
    }
}


int
main()
{
  initlog();
  check<3, 3>();
}
