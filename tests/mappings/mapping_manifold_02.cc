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


// Check that the quadrature points obtained through transform_unit_to_real
// are the same of those given by FEValues


#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>

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
  deallog << "dim=" << dim << ", spacedim=" << spacedim << std::endl;

  Triangulation<dim, spacedim> triangulation;

  Point<spacedim> center;
  for (unsigned int i = 0; i < spacedim; ++i)
    center[i] = 5 + i;

  const double inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(triangulation, center, inner_radius, outer_radius);


  static const PolarManifold<dim> manifold(center);

  triangulation.set_all_manifold_ids(0);
  triangulation.set_manifold(0, manifold);

  const QGauss<spacedim> quad(3);

  FE_Q<dim, spacedim> fe(1);

  DoFHandler<dim, spacedim> dof(triangulation);
  dof.distribute_dofs(fe);

  MappingManifold<dim, spacedim> mapping;

  FEValues<dim, spacedim> fe_values(mapping,
                                    fe,
                                    quad,
                                    update_quadrature_points);

  for (typename Triangulation<dim, spacedim>::active_cell_iterator cell =
         triangulation.begin_active();
       cell != triangulation.end();
       ++cell)
    {
      fe_values.reinit(cell);
      std::vector<Point<spacedim>> fev_qp = fe_values.get_quadrature_points();
      for (unsigned int q = 0; q < fev_qp.size(); ++q)
        {
          const Point<spacedim> pq =
            mapping.transform_unit_to_real_cell(cell, quad.point(q));

          if (pq.distance(fev_qp[q]) > 1e-10)
            {
              deallog << "Expected: " << pq << ", got: " << fev_qp[q]
                      << std::endl;
            }
        }
    }
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test<2, 2>();
  test<3, 3>();

  return 0;
}
