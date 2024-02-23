// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that MappingManifold and MappingQ1 are the same thing on a
// FlatManifold. Test on the quadrature points.

#include <deal.II/base/utilities.h>

#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  deallog << "dim=" << dim << ", spacedim=" << spacedim << std::endl;

  Triangulation<dim, spacedim> triangulation;

  GridGenerator::hyper_cube(triangulation, 2.0, 3.0);


  const QGauss<dim> quadrature(5);

  std::vector<Point<dim>> q_points = quadrature.get_points();

  MappingManifold<dim, spacedim> map_manifold;
  MappingQ<dim, spacedim>        map_q1(1);

  typename Triangulation<dim, spacedim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
  for (; cell != endc; ++cell)
    {
      for (unsigned int i = 0; i < q_points.size(); ++i)
        {
          const Point<spacedim> pq =
            map_manifold.transform_unit_to_real_cell(cell, q_points[i]);
          const Point<spacedim> pq1 =
            map_q1.transform_unit_to_real_cell(cell, q_points[i]);
          if (pq.distance(pq1) > 1e-10)
            {
              deallog << "Expected: " << pq << ", got: " << pq1 << std::endl;
            }
        }
    }
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test<1, 1>();
  test<2, 2>();
  test<3, 3>();

  test<1, 2>();
  test<1, 3>();
  test<2, 3>();

  return 0;
}
