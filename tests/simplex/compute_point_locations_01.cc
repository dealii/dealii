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

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include "../tests.h"

// Create a simplex mesh in the unit cube. Check points distribute to each cell
// via GridTools::compute_point_locations()


template <int dim>
void
test_in_unit_cube(const std::vector<Point<dim>> &points)
{
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);

  GridTools::Cache<dim> cache(tria);
  const auto            point_locations =
    GridTools::compute_point_locations(cache, points);

  const auto cells   = std::get<0>(point_locations);
  const auto qpoints = std::get<1>(point_locations);
  const auto indices = std::get<2>(point_locations);

  deallog << "Print out results in " << dim << "D: " << std::endl;
  for (unsigned int i = 0; i < cells.size(); ++i)
    {
      deallog << "At cell " << i << ':' << std::endl;
      for (unsigned int j = 0; j < qpoints[i].size(); ++j)
        {
          deallog << "    qpoints " << j << ": " << std::endl;
          deallog << "        reference position  : (" << qpoints[i][j] << ')'
                  << std::endl;
          deallog << "        physical position   : (" << points[indices[i][j]]
                  << ')' << std::endl;
          deallog << "        FE mapping position : ("
                  << cache.get_mapping().transform_unit_to_real_cell(
                       cells[i], qpoints[i][j])
                  << ')' << std::endl;
        }
    }
  deallog << std::endl;
}

int
main()
{
  initlog();

  {
    unsigned int          n = 3;
    std::vector<Point<2>> test_points((n + 1) * (n + 1));
    unsigned int          global_index = 0;
    for (unsigned int i = 0; i < n + 1; ++i)
      for (unsigned int j = 0; j < n + 1; ++j)
        {
          test_points[global_index][0] = (1.0 / n) * i;
          test_points[global_index][1] = (1.0 / n) * j;
          ++global_index;
        }
    test_in_unit_cube(test_points);
  }

  {
    unsigned int          n = 3;
    std::vector<Point<3>> test_points((n + 1) * (n + 1) * (n + 1));
    unsigned int          global_index = 0;
    for (unsigned int i = 0; i < n + 1; ++i)
      for (unsigned int j = 0; j < n + 1; ++j)
        for (unsigned int k = 0; k < n + 1; ++k)
          {
            test_points[global_index][0] = (1.0 / n) * i;
            test_points[global_index][1] = (1.0 / n) * j;
            test_points[global_index][2] = (1.0 / n) * k;
            ++global_index;
          }
    test_in_unit_cube(test_points);
  }
}
