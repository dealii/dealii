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

// Try to compute the area of a circle/sphere using JxW values.

#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
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
  std::ostream &out = deallog.get_file_stream();

  out << "# dim=" << dim << ", spacedim=" << spacedim << std::endl;

  Triangulation<dim, spacedim> triangulation;

  Point<spacedim> center;
  center[0] = 1.5;
  center[1] = 2.5;

  double radius = 1.0;

  static const PolarManifold<dim, spacedim> manifold(center);
  GridGenerator::hyper_ball(triangulation, center, radius);

  triangulation.set_all_manifold_ids_on_boundary(0);
  triangulation.set_manifold(0, manifold);

  for (unsigned int cycle = 0; cycle < 4; ++cycle)
    {
      MappingManifold<dim, spacedim> map_manifold;
      FE_Q<dim, spacedim>            fe(1);
      const QGauss<dim - 1>          quad(3);

      FEFaceValues<dim, spacedim> fe_v(map_manifold,
                                       fe,
                                       quad,
                                       update_JxW_values);
      double                      area = 0;

      for (typename Triangulation<dim, spacedim>::active_cell_iterator cell =
             triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        for (const unsigned int f : GeometryInfo<dim>::face_indices())
          if (cell->face(f)->at_boundary())
            {
              fe_v.reinit(cell, f);
              for (unsigned int i = 0; i < quad.size(); ++i)
                area += fe_v.JxW(i);
            }
      deallog << "Cycle	      : " << cycle << std::endl;
      deallog << "Surface Area  : " << area << std::endl;
      deallog << "Error         : " << (area - (dim - 1) * 2 * numbers::PI)
              << std::endl;

      triangulation.refine_global(1);
    }
}


int
main()
{
  initlog();

  test<2, 2>();

  test<3, 3>();

  return 0;
}
