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

// Compute intersection of a segment embedded in 2D and a deal.II cell in 2D,
// and return a vector of arrays where you can build Quadrature rules. Then
// check that the sum of weights give the correct area for each region.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/cgal/intersections.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

using namespace CGALWrappers;
static QGaussSimplex<1> qgauss(1); // use a degree 1 QGaussSimplex<2>

void
test_inside_intersection(Triangulation<2> &tria0, Triangulation<1, 2> &tria1)
{
  GridGenerator::hyper_cube(tria0, -1., 1.);
  GridGenerator::hyper_cube(tria1, -0.45, 0.45);
  const double expected_area = 0.9;
  const auto   cell0         = tria0.begin_active();
  const auto   cell1         = tria1.begin_active();
  // cell1->vertex(1)+=Point<2>(1.0,1.);

  const auto vec_of_arrays = CGALWrappers::compute_intersection_of_cells(
    cell0, cell1, MappingQ1<2>(), MappingQ1<1, 2>());

  const auto   quad = qgauss.mapped_quadrature(vec_of_arrays);
  const double sum =
    std::accumulate(quad.get_weights().begin(), quad.get_weights().end(), 0.);

  assert(std::abs(expected_area - sum) < 1e-15);
  deallog << "OK" << std::endl;
}



void
test_intersection(Triangulation<2> &tria0, Triangulation<1, 2> &tria1)
{
  GridGenerator::hyper_cube(tria0, -1., 1.);
  GridGenerator::hyper_cube(tria1, .5, 1.45);
  const double expected_area = .5;
  const auto   cell0         = tria0.begin_active();
  const auto   cell1         = tria1.begin_active();

  const auto vec_of_arrays = CGALWrappers::compute_intersection_of_cells(
    cell0, cell1, MappingQ1<2>(), MappingQ1<1, 2>());

  const auto   quad = qgauss.mapped_quadrature(vec_of_arrays);
  const double sum =
    std::accumulate(quad.get_weights().begin(), quad.get_weights().end(), 0.);
  assert(std::abs(expected_area - sum) < 1e-15);
  deallog << "OK" << std::endl;
}



void
test_failing_intersection(Triangulation<2> &tria0, Triangulation<1, 2> &tria1)
{
  GridGenerator::hyper_cube(tria0, -1., 1.);
  GridGenerator::hyper_cube(tria1, 1.5, 2.);
  const double expected_area = 0.;
  const auto   cell0         = tria0.begin_active();
  const auto   cell1         = tria1.begin_active();

  const auto vec_of_arrays = CGALWrappers::compute_intersection_of_cells(
    cell0, cell1, MappingQ1<2>(), MappingQ1<1, 2>());

  const auto   quad = qgauss.mapped_quadrature(vec_of_arrays);
  const double sum =
    std::accumulate(quad.get_weights().begin(), quad.get_weights().end(), 0.);
  assert(std::abs(expected_area - sum) < 1e-15);
  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  Triangulation<2>    tria0; // ambient space
  Triangulation<1, 2> tria1; // immersed grid
  test_inside_intersection(tria0, tria1);

  tria0.clear();
  tria1.clear();

  test_intersection(tria0, tria1);

  tria0.clear();
  tria1.clear();

  test_failing_intersection(tria0, tria1);
}
