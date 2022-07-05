// ---------------------------------------------------------------------

// Copyright (C) 2022 by the deal.II authors

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

// Compute intersection of two 3D cells, and return a Quadrature rule over it.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/cgal/intersections.h>

#include "../tests.h"

static QGaussSimplex<3> qgauss(1);
void
test_intersection(Triangulation<3> &tria0, Triangulation<3> &tria1)
{
  GridGenerator::hyper_cube(tria0, -1., 1.);
  GridGenerator::hyper_cube(tria1, .5, 1.5);
  static constexpr double expected_measure = 0.5 * 0.5 * 0.5;
  const auto              cell0            = tria0.begin_active();
  const auto              cell1            = tria1.begin_active();

  const auto vec_of_arrays =
    CGALWrappers::compute_intersection_of_cells<3, 3, 3>(cell0,
                                                         cell1,
                                                         MappingQ1<3>(),
                                                         MappingQ1<3>());

  const auto   quad = qgauss.mapped_quadrature(vec_of_arrays);
  const double sum =
    std::accumulate(quad.get_weights().begin(), quad.get_weights().end(), 0.);
  assert(std::abs(sum - expected_measure) < 1e-15);
  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  Triangulation<3> tria0;
  Triangulation<3> tria1;

  test_intersection(tria0, tria1);
}
