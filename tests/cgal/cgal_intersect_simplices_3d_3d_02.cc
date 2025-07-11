// ------------------------------------------------------------------------

// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors

// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Compute intersection of two 3D cells. This additional test is added because
// intersections are not found with inexact kernels.

#include <deal.II/cgal/intersections.h>
#include <deal.II/cgal/utilities.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

void
test_intersection()
{
  MappingQ1<3> mapping;

  Triangulation<3> tria0;
  GridGenerator::hyper_cube(tria0);
  const auto cell0 = tria0.begin_active();
  {
    cell0->vertex(0) =
      Point<3>(2.600000000000e-01, 0.000000000000e+00, 0.000000000000e+00);
    cell0->vertex(1) =
      Point<3>(2.600000000000e-01, 0.000000000000e+00, 6.955875000000e-03);
    cell0->vertex(2) =
      Point<3>(2.600000000000e-01, 6.955875000000e-03, 0.000000000000e+00);
    cell0->vertex(3) =
      Point<3>(2.600000000000e-01, 6.158125000000e-03, 6.158125000000e-03);
    cell0->vertex(4) =
      Point<3>(2.550000000000e-01, 0.000000000000e+00, 0.000000000000e+00);
    cell0->vertex(5) =
      Point<3>(2.550000000000e-01, 0.000000000000e+00, 6.955875000000e-03);
    cell0->vertex(6) =
      Point<3>(2.550000000000e-01, 6.955875000000e-03, 0.000000000000e+00);
    cell0->vertex(7) =
      Point<3>(2.550000000000e-01, 6.158125000000e-03, 6.158125000000e-03);
  }

  Triangulation<3> tria1;
  GridGenerator::hyper_cube(tria1);
  const auto cell1 = tria1.begin_active();
  {
    cell1->vertex(0) =
      Point<3>(2.600000000000e-01, 0.000000000000e+00, 0.000000000000e+00);
    cell1->vertex(1) =
      Point<3>(2.600000000000e-01, 0.000000000000e+00, 1.391175000000e-02);
    cell1->vertex(2) =
      Point<3>(2.600000000000e-01, 1.391175000000e-02, 0.000000000000e+00);
    cell1->vertex(3) =
      Point<3>(2.600000000000e-01, 1.072075000000e-02, 1.072075000000e-02);
    cell1->vertex(4) =
      Point<3>(2.500000000000e-01, 0.000000000000e+00, 0.000000000000e+00);
    cell1->vertex(5) =
      Point<3>(2.500000000000e-01, 0.000000000000e+00, 1.391175000000e-02);
    cell1->vertex(6) =
      Point<3>(2.500000000000e-01, 1.391175000000e-02, 0.000000000000e+00);
    cell1->vertex(7) =
      Point<3>(2.500000000000e-01, 1.072075000000e-02, 1.072075000000e-02);
  }

  const auto &vec_of_arrays = CGALWrappers::compute_intersection_of_cells(
    cell1, cell0, mapping, mapping, 1e-9);

  assert(vec_of_arrays.size() == 18);

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  test_intersection();

  return 0;
}
