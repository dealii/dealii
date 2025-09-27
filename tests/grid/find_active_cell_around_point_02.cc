// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/*
 * Test cached version of find active cell around point
 */
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test(unsigned int n_ref, unsigned int n_points)
{
  deallog << "Testing " << dim << ", " << spacedim << std::endl;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_ref);
  DoFHandler<dim, spacedim> dof_handler(tria);

  std::vector<Point<spacedim>> points;

  deallog << "Points in study: " << n_points << std::endl;
  for (std::size_t i = 0; i < n_points; ++i)
    points.push_back(random_point<spacedim>());

  auto &mapping = StaticMappingQ1<dim, spacedim>::mapping;

  GridTools::Cache<dim, spacedim> cache(tria, mapping);

  auto cell = tria.begin_active();
  for (auto &p : points)
    {
      auto c_and_p = GridTools::find_active_cell_around_point(cache, p);
      auto p2 =
        mapping.transform_unit_to_real_cell(c_and_p.first, c_and_p.second);
      if (p2.distance(p) > 1e-10)
        deallog << "NOT OK!" << p << ", " << p2 << ", " << c_and_p.first
                << std::endl;
      cell = c_and_p.first;
    }
  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  test<1, 1>(4, 10);
  test<2, 2>(3, 20);
  test<3, 3>(2, 30);
}
