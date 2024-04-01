// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check MappingCartesian::transform_points_real_to_unit_cell, using a case
// otherwise similar to mapping_points_real_to_unit.cc

#include <deal.II/base/utilities.h>

#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  Point<dim>         lower, upper;
  for (unsigned int d = 0; d < dim; ++d)
    lower[d] = -0.3 + d * 0.07;
  for (unsigned int d = 0; d < dim; ++d)
    upper[d] = 0.45 - d * 0.03;

  GridGenerator::hyper_rectangle(triangulation, lower, upper);

  deallog << "dim=" << dim << std::endl;
  deallog << "MappingCartesian: ";
  const unsigned int      n_points = 5;
  std::vector<Point<dim>> unit_points(Utilities::pow(n_points, dim));

  switch (dim)
    {
      case 1:
        for (unsigned int x = 0; x < n_points; ++x)
          unit_points[x][0] = static_cast<double>(x) / n_points;
        break;

      case 2:
        for (unsigned int x = 0; x < n_points; ++x)
          for (unsigned int y = 0; y < n_points; ++y)
            {
              unit_points[y * n_points + x][0] =
                static_cast<double>(x) / n_points;
              unit_points[y * n_points + x][1] =
                static_cast<double>(y) / n_points;
            }
        break;

      case 3:
        for (unsigned int x = 0; x < n_points; ++x)
          for (unsigned int y = 0; y < n_points; ++y)
            for (unsigned int z = 0; z < n_points; ++z)
              {
                unit_points[z * n_points * n_points + y * n_points + x][0] =
                  static_cast<double>(x) / n_points;
                unit_points[z * n_points * n_points + y * n_points + x][1] =
                  static_cast<double>(y) / n_points;
                unit_points[z * n_points * n_points + y * n_points + x][2] =
                  static_cast<double>(z) / n_points;
              }
        break;
    }

  MappingCartesian<dim>                             mapping;
  typename Triangulation<dim>::active_cell_iterator cell =
    triangulation.begin_active();
  std::vector<Point<dim>> real_points(unit_points.size());
  for (unsigned int i = 0; i < unit_points.size(); ++i)
    real_points[i] = mapping.transform_unit_to_real_cell(cell, unit_points[i]);
  std::vector<Point<dim>> new_points(unit_points.size());
  mapping.transform_points_real_to_unit_cell(cell, real_points, new_points);
  for (unsigned int i = 0; i < unit_points.size(); ++i)
    {
      // for each of the points, verify that applying the forward map and
      // then pull back get the same point again
      AssertThrow(unit_points[i].distance(new_points[i]) < 1e-10,
                  ExcInternalError());
    }
  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
