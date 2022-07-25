// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2022 by the deal.II authors
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



// check GridTools::cell_measure()

#include <deal.II/grid/grid_tools.h>

#include "../tests.h"


template <int dim>
void
test()
{}

template <>
void
test<1>()
{
  std::vector<Point<1>> points;
  points.emplace_back(0.);
  points.emplace_back(3.);

  unsigned int indices1[GeometryInfo<1>::vertices_per_cell] = {0, 1};
  std::vector<unsigned int> indices2                        = {0, 1};

  deallog << "1d: " << GridTools::cell_measure(points, indices1) << ' '
          << GridTools::cell_measure(points, indices2) << std::endl;
}

template <>
void
test<2>()
{
  std::vector<Point<2>> points;
  points.emplace_back(0., 0.);
  points.emplace_back(1., 0.);
  points.emplace_back(0., 2.);
  points.emplace_back(1., 2.);

  unsigned int indices1[GeometryInfo<2>::vertices_per_cell] = {0, 1, 2, 3};
  std::vector<unsigned int> indices2                        = {0, 1, 2, 3};

  deallog << "2d: " << GridTools::cell_measure(points, indices1) << ' '
          << GridTools::cell_measure(points, indices2) << std::endl;
}

template <>
void
test<3>()
{
  std::vector<Point<3>> points;
  points.emplace_back(0., 0., 0.);
  points.emplace_back(1., 0., 0.);
  points.emplace_back(0., 2., 0.);
  points.emplace_back(1., 2., 0.);
  points.emplace_back(0., 0., 3.);
  points.emplace_back(1., 0., 3.);
  points.emplace_back(0., 2., 3.);
  points.emplace_back(1., 2., 3.);

  unsigned int indices1[GeometryInfo<3>::vertices_per_cell] = {
    0, 1, 2, 3, 4, 5, 6, 7};
  std::vector<unsigned int> indices2 = {0, 1, 2, 3, 4, 5, 6, 7};

  deallog << "3d: " << GridTools::cell_measure(points, indices1) << ' '
          << GridTools::cell_measure(points, indices2) << std::endl;
}

int
main()
{
  initlog();
  test<1>();
  test<2>();
  test<3>();

  return 0;
}
