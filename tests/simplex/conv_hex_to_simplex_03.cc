/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 - 2021 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 * Like conv_hex_to_simplex, but also refines the grid once.
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
