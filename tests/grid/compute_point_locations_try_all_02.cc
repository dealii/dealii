// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2023 by the deal.II authors
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

// Test GridTools::compute_point_locations_try_all for the case where
// some points only lie outside a cell by some roundoff, and another case
// where they are farther away

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
check_found(const Triangulation<dim>      &tria,
            const GridTools::Cache<dim>   &cache,
            const std::vector<Point<dim>> &points)
{
  auto cell_qpoint_map =
    GridTools::compute_point_locations_try_all(cache, points);
  const auto &cells   = std::get<0>(cell_qpoint_map);
  const auto &qpoints = std::get<1>(cell_qpoint_map);
  const auto &indices = std::get<2>(cell_qpoint_map);
  const auto &other_p = std::get<3>(cell_qpoint_map);
  size_t      n_cells = cells.size();

  deallog << "Points found in " << n_cells << " cells" << std::endl;

  for (unsigned int i = 0; i < n_cells; ++i)
    {
      deallog << "On cell " << cells[i]->id() << " found:" << std::endl;
      unsigned int j = 0;
      for (const Point<dim> &p : qpoints[i])
        deallog << "real " << points[indices[i][j++]] << " unit " << p
                << std::endl;
    }
  deallog << "Points not found: ";
  for (const auto i : other_p)
    deallog << points[i] << "   ";
  deallog << std::endl << std::endl;
}


int
main()
{
  initlog();

  deallog << std::setprecision(18);

  const int          dim = 2;
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -0.3, 0.7);
  triangulation.refine_global(2);

  MappingQ<dim>         mapping(1);
  GridTools::Cache<dim> cache(triangulation, mapping);

  // Point is outside by 1e-16
  Point<dim> p0(-0.299999999999999989, -0.247168783648703205),
    p1(-0.300000000000000044, -0.102831216351296786);

  check_found(triangulation, cache, {{p0, p1}});

  // Shift point by 1e-8
  p0[0] -= 1e-8;

  check_found(triangulation, cache, {{p0, p1}});

  // Shift point by more
  p0[0] += 0.5;
  p1[0] += 1e-15;

  check_found(triangulation, cache, {{p0, p1}});

  // Shift outside
  p0[0] += 1.;

  check_found(triangulation, cache, {{p0, p1}});
}
