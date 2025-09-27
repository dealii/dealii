// ------------------------------------------------------------------------

// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors

// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Compute intersection of two 3D cells, and return a Quadrature rule over it.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/cgal/intersections.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

void
test_intersection(Triangulation<3> &tria0, Triangulation<3> &tria1)
{
  GridGenerator::hyper_cube(tria0, -1., 1.);
  GridGenerator::hyper_cube(tria1, .5, 1.5);
  static constexpr double expected_measure = 0.5 * 0.5 * 0.5;
  const auto              cell0            = tria0.begin_active();
  const auto              cell1            = tria1.begin_active();

  const auto vec_of_arrays = CGALWrappers::compute_intersection_of_cells(
    cell0, cell1, MappingQ1<3>(), MappingQ1<3>());

  const auto   quad = QGaussSimplex<3>(1).mapped_quadrature(vec_of_arrays);
  const double sum =
    std::accumulate(quad.get_weights().begin(), quad.get_weights().end(), 0.);
  assert(std::abs(sum - expected_measure) < 1e-15);
  deallog << "OK" << std::endl;
}



void
test_failing_intersection(Triangulation<3> &tria0, Triangulation<3> &tria1)
{
  GridGenerator::hyper_cube(tria0, -1., 1.);
  GridGenerator::hyper_cube(tria1, 2.5, 3.);
  static constexpr double expected_measure = 0.;
  const auto              cell0            = tria0.begin_active();
  const auto              cell1            = tria1.begin_active();

  const auto vec_of_arrays = CGALWrappers::compute_intersection_of_cells(
    cell0, cell1, MappingQ1<3>(), MappingQ1<3>());

  const auto   quad = QGaussSimplex<3>(1).mapped_quadrature(vec_of_arrays);
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

  tria0.clear();
  tria1.clear();

  test_failing_intersection(tria0, tria1);
}
