// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Step-01 on a simplex mesh. Following incompatible modifications had to be
// made:
//  - Replace or convert triangulation to a simplex based triangulation.
//  - Copy manifolds after conversion to simplex based triangulation in order
//    to execute refinement.


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <cmath>
#include <fstream>
#include <iostream>

#include "../tests.h"


void
first_grid()
{
  Triangulation<2> triangulation;
  GridGenerator::subdivided_hyper_cube_with_simplices(triangulation, 1);
  triangulation.refine_global(2);

  GridOut grid_out;
  grid_out.write_svg(triangulation, deallog.get_file_stream());
  deallog << "Grid written to grid-1.svg" << std::endl;
}



void
second_grid()
{
  Triangulation<2> tria_temp, triangulation;

  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(tria_temp, center, inner_radius, outer_radius, 10);

  // convert mesh to a simplex mesh and copy manually the attached manifolds
  GridGenerator::convert_hypercube_to_simplex_mesh(tria_temp, triangulation);
  for (const auto i : tria_temp.get_manifold_ids())
    if (i != numbers::flat_manifold_id)
      triangulation.set_manifold(i, tria_temp.get_manifold(i));

  for (unsigned int step = 0; step < 2; ++step)
    {
      for (auto &cell : triangulation.active_cell_iterators())
        {
          for (const auto v : cell->vertex_indices())
            {
              const double distance_from_center =
                center.distance(cell->vertex(v));

              if (std::fabs(distance_from_center - inner_radius) <=
                  1e-6 * inner_radius)
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }

      triangulation.execute_coarsening_and_refinement();
    }

  GridOut grid_out;
#if false
  std::ofstream out("mesh.svg");
  grid_out.write_svg(triangulation, out);
#else
  grid_out.write_svg(triangulation, deallog.get_file_stream());
#endif

  deallog << "Grid written to grid-2.svg" << std::endl;
}



int
main()
{
  initlog();

  first_grid();
  second_grid();
}
