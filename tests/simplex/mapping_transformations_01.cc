// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// Create a simplex mesh, and transform back and forth a point in each cell

#include <deal.II/base/patterns.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <fstream>

#include "../tests.h"

using namespace dealii;

void make_grid(Triangulation<2> &triangulation)
{
  Triangulation<2> triangulation_temp;

  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation_temp, center, inner_radius, outer_radius, 5);

  GridGenerator::convert_hypercube_to_simplex_mesh(triangulation_temp,
                                                   triangulation);
  triangulation.set_manifold(0, FlatManifold<2>());
}

int
main()
{
  initlog();

  Triangulation<2> triangulation;
  make_grid(triangulation);

  MappingFE<2> mapping(FE_SimplexP<2>(1));

  unsigned int n_points = 1;

  std::vector<Point<2>> points;
  do
    {
      auto p   = random_point<2>();
      auto sum = p[0] + p[1];
      if (sum <= 1)
        points.push_back(p);
    }
  while (points.size() < n_points);

  deallog << "Transforming " << n_points << " from reference to real."
          << std::endl
          << "Points: " << Patterns::Tools::to_string(points) << std::endl;

  for (const auto &cell : triangulation.active_cell_iterators())
    for (const auto &p : points)
      {
        const auto real_p = mapping.transform_unit_to_real_cell(cell, p);
        const auto pull_back_p =
          mapping.transform_real_to_unit_cell(cell, real_p);
        deallog << "F(" << p << ")=" << real_p << " --- "
                << "F^-1(" << real_p << ")=" << pull_back_p << std::endl;
        if (p.distance(pull_back_p) > 1e-10)
          deallog << "ERROR: " << p.distance(pull_back_p) << std::endl;
      }
}
