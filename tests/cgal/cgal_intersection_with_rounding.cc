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

// Compute intersection of simplices in 2D, and return a vector of arrays of
// intersections. If elements are not intersecting due to rounding errors we
// specify a parameter where to truncate the vertices of the elements.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/cgal/intersections.h>

#include "../tests.h"

using namespace CGALWrappers;

void
test_truncated_intersection(Triangulation<2, 2> &tria0,
                            Triangulation<1, 2> &tria1,
                            const double         eps0,
                            const double         eps1)
{
  GridGenerator::hyper_cube(tria0);
  GridGenerator::hyper_cube(tria1);

  // cells
  //|------||
  //|   0  || 1
  //|------||

  const auto cell0 = tria0.begin_active();
  cell0->vertex(0) = Point<2>(0.0, 0.0);
  cell0->vertex(1) = Point<2>(1.0 - eps0, 0.0);
  cell0->vertex(2) = Point<2>(0.0, 1.0);
  cell0->vertex(3) = Point<2>(1.0 - eps0, 1.0);

  const auto cell1 = tria1.begin_active();
  cell1->vertex(0) = Point<2>(1.0 + eps1, 0.0);
  cell1->vertex(1) = Point<2>(1.0 + eps1, 1.0);

  {
    const auto vec_of_arrays =
      CGALWrappers::compute_intersection_of_cells<2, 1, 2>(
        cell0, cell1, MappingQ1<2, 2>(), MappingQ1<1, 2>(), 1e-9);
    deallog << "number of found intersections: " << vec_of_arrays.size()
            << std::endl;
    deallog << "OK" << std::endl;
  }
  {
    const auto vec_of_arrays =
      CGALWrappers::compute_intersection_of_cells<2, 1, 2>(
        cell0, cell1, MappingQ1<2, 2>(), MappingQ1<1, 2>(), 1e-9, 1e-11);
    deallog << "number of found intersections: " << vec_of_arrays.size()
            << std::endl;
    deallog << "OK" << std::endl;
  }
}


int
main()
{
  initlog();
  Triangulation<2, 2> tria0;
  Triangulation<1, 2> tria1;

  // cell0_x0=1.0, cell1_x0=1.0
  test_truncated_intersection(tria0, tria1, 0.0, 0.0);

  tria0.clear();
  tria1.clear();

  // cell0_x0=1.0, cell1_x0=1.0 + 1e-12
  test_truncated_intersection(tria0, tria1, 0.0, 1e-12);

  tria0.clear();
  tria1.clear();

  // cell0_x0=1.0 - 1e-12, cell1_x0=1.0
  test_truncated_intersection(tria0, tria1, 1e-12, 0.0);
}
