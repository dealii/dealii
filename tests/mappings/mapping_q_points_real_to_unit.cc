// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check MappingQ::transform_real_to_unit_points on a set of
// challenging points, especially with vectorization because we have nearby
// points that succeed and others in the regime of negative Jacobian
// determinants, respectively

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test(const unsigned int degree)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_shell(tria, Point<dim>(), 0.5, 1., 2 * dim);

  std::vector<Point<dim>> real_points;
  Point<dim>              p;
  for (unsigned int d = 0; d < dim; ++d)
    p[d] = std::sqrt(1. / dim) + 0.01 * d;
  for (unsigned i = 1; i < 21; ++i)
    real_points.push_back(
      (i % 2 ? 0.05 * (static_cast<double>(i) - 10.) : 0.05 * i) * p);
  std::vector<Point<dim>> unit_points(real_points.size());
  MappingQ<dim>           mapping(degree);
  mapping.transform_points_real_to_unit_cell(tria.begin(),
                                             real_points,
                                             unit_points);

  deallog << "Transform on cell with center: " << tria.begin()->center(true)
          << " with mapping degree " << degree << std::endl;
  for (unsigned int i = 0; i < real_points.size(); ++i)
    deallog << "Transform " << real_points[i] << " gives " << unit_points[i]
            << std::endl;
  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  test<2>(5);
}
