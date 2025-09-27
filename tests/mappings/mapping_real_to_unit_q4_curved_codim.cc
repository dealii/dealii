// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// on a somewhat deformed cube, verify that if we push forward a bunch
// of points from the reference to the real cell and then call
// Mapping::transform_unit_to_real_cell that we get the same point as
// we had in the beginning.
//
// like in the _q4_straight test, we use a Q4 mapping but this time we
// actually curve one boundary of the cell which ensures that the
// mapping is really higher order than just Q1

#include <deal.II/base/utilities.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test_real_to_unit_cell()
{
  deallog << "dim=" << dim << ", spacedim=" << spacedim << std::endl;

  // define a boundary that fits the
  // the vertices of the hyper cube
  // we're going to create below
  SphericalManifold<dim, spacedim> boundary;

  Triangulation<dim, spacedim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);

  // set the boundary indicator for
  // one face of the single cell
  triangulation.set_manifold(1, boundary);
  triangulation.begin_active()->face(0)->set_boundary_id(1);

  const unsigned int      n_points = 5;
  std::vector<Point<dim>> unit_points(Utilities::fixed_power<dim>(n_points));

  switch (dim)
    {
      case 1:
        for (unsigned int x = 0; x < n_points; ++x)
          unit_points[x][0] = double(x) / double(n_points);
        break;

      case 2:
        for (unsigned int x = 0; x < n_points; ++x)
          for (unsigned int y = 0; y < n_points; ++y)
            {
              unit_points[y * n_points + x][0] = double(x) / double(n_points);
              unit_points[y * n_points + x][1] = double(y) / double(n_points);
            }
        break;

      case 3:
        for (unsigned int x = 0; x < n_points; ++x)
          for (unsigned int y = 0; y < n_points; ++y)
            for (unsigned int z = 0; z < n_points; ++z)
              {
                unit_points[z * n_points + y * n_points + x][0] =
                  double(x) / double(n_points);
                unit_points[z * n_points + y * n_points + x][1] =
                  double(y) / double(n_points);
                unit_points[z * n_points + y * n_points + x][2] =
                  double(z) / double(n_points);
              }
        break;
    }


  MappingQ<dim, spacedim> map(4);

  // work with this cell (unlike the
  // _q1 test where we move vertices)
  typename Triangulation<dim, spacedim>::active_cell_iterator cell =
    triangulation.begin_active();
  for (unsigned int i = 0; i < unit_points.size(); ++i)
    {
      // for each of the points,
      // verify that if we apply
      // the forward map and then
      // pull back that we get
      // the same point again
      const Point<spacedim> p =
        map.transform_unit_to_real_cell(cell, unit_points[i]);
      const Point<dim> p_unit = map.transform_real_to_unit_cell(cell, p);
      AssertThrow(unit_points[i].distance(p_unit) < 1e-10, ExcInternalError());
    }
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test_real_to_unit_cell<2, 3>();
}
