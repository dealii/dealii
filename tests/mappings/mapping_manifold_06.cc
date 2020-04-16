// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

// Check MappingManifold on faces. Create a sphere and make sure all
// points returned by FEValues.get_quadrature_points are actually on
// the circle.

#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/mapping_q_generic.h>

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

  MappingManifold<dim, spacedim> map_manifold;
  FE_Q<dim, spacedim>            fe(1);
  const QGauss<dim - 1>          quad(5);

  FEFaceValues<dim, spacedim> fe_v(
    map_manifold, fe, quad, update_quadrature_points | update_normal_vectors);

  out << "set size ratio -1" << std::endl
      << "set terminal aqua " << dim + spacedim << std::endl
      << (spacedim == 3 ? "s" : "") << "plot '-' w l" << std::endl;

  for (typename Triangulation<dim, spacedim>::active_cell_iterator cell =
         triangulation.begin_active();
       cell != triangulation.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      {
        const double xc = cell->face(f)->center()[0] - 1.5;
        const double yc = cell->face(f)->center()[1] - 2.5;
        const double rc = xc * xc + yc * yc;

        if (cell->face(f)->at_boundary() && rc > 0.1)
          {
            fe_v.reinit(cell, f);
            const std::vector<Point<spacedim>> &qps =
              fe_v.get_quadrature_points();
            const std::vector<Tensor<1, spacedim>> &nps =
              fe_v.get_normal_vectors();
            for (unsigned int i = 0; i < qps.size(); ++i)
              {
                out << qps[i] << std::endl;
                if (std::abs(qps[i].distance(center) - radius) > 1e-10)
                  out << "# Error! This should be on the sphere, but it's not!"
                      << std::endl;
              }
            out << std::endl;
          }
      }

  out << "e" << std::endl
      << std::endl
      << std::endl
      << "set terminal aqua " << 2 * (dim + spacedim) + 1 << std::endl
      << (spacedim == 3 ? "s" : "") << "plot '-' with vectors" << std::endl;

  for (typename Triangulation<dim, spacedim>::active_cell_iterator cell =
         triangulation.begin_active();
       cell != triangulation.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      {
        const double xc = cell->face(f)->center()[0] - 1.5;
        const double yc = cell->face(f)->center()[1] - 2.5;
        const double rc = xc * xc + yc * yc;

        if (cell->face(f)->at_boundary() && rc > 0.1)
          {
            fe_v.reinit(cell, f);
            const std::vector<Point<spacedim>> &qps =
              fe_v.get_quadrature_points();
            const std::vector<Tensor<1, spacedim>> &nps =
              fe_v.get_normal_vectors();
            for (unsigned int i = 0; i < qps.size(); ++i)
              {
                out << qps[i] << " " << nps[i] << std::endl;
                if (std::abs((center + nps[i] - qps[i]).norm()) > 1e-10)
                  out << "# Error! The normal vector should be radial! "
                      << std::endl
                      << "# Got " << center + nps[i] << " instead of " << qps[i]
                      << std::endl;
              }
            out << std::endl;
          }
      }
  out << "e" << std::endl << std::endl;
}


int
main()
{
  initlog();

  test<2, 2>();

  test<3, 3>();

  return 0;
}
