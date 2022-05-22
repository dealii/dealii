// Copyright (C) 2020 - 2022 by the deal.II authors
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

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include "../tests.h"

// Create a simplex mesh in the unit cube. Check points distribute to each cell
// via GridTools::compute_point_locations

using namespace dealii;

template <int dim>
void
test_in_unit_cube(const std::vector<Point<dim>> &points)
{
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);

  MappingFE<dim> mapping(FE_SimplexP<dim>(1));

  const auto tria_cache =
    std::make_unique<GridTools::Cache<dim>>(tria, mapping);

  const auto point_locations =
    GridTools::compute_point_locations(*tria_cache, points);

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
                  << mapping.transform_unit_to_real_cell(cells[i],
                                                         qpoints[i][j])
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
