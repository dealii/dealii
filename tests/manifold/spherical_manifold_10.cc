// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Compare SphericalManifold::normal_vector with the result of
// FEFaceValues::normal_vector. Those should give very similar results, but an
// old implementation in the manifold gave wrong points

#include <deal.II/base/point.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(8);

  constexpr unsigned int dim = 3;
  SphericalManifold<3>   spherical;

  Triangulation<dim> tria;
  GridGenerator::hyper_shell(tria, Point<dim>(), 0.5, 1., 96, true);
  tria.set_all_manifold_ids(0);
  tria.set_manifold(0, spherical);

  MappingQ<dim>          mapping(4);
  QGaussLobatto<dim - 1> quadrature(4);

  FE_Nothing<dim>   dummy;
  FEFaceValues<dim> fe_values(mapping,
                              dummy,
                              quadrature,
                              update_normal_vectors | update_quadrature_points);

  for (auto &cell : tria.active_cell_iterators())
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      {
        fe_values.reinit(cell, f);

        // all points should coincide to an accuracy of at least 1e-8 (note
        // that we use a 4-th degree mapping, so its accuracy on the 96 cell
        // version of the sphere should be enough)
        const double tolerance = 1e-8;
        for (unsigned int q = 0; q < quadrature.size(); ++q)
          {
            const Tensor<1, dim> normal_manifold =
              spherical.normal_vector(cell->face(f),
                                      fe_values.quadrature_point(q));
            const Tensor<1, dim> normal_feval = fe_values.normal_vector(q);
            if (std::abs(1.0 - std::abs(normal_manifold * normal_feval)) >
                tolerance)
              deallog << "Error in point " << fe_values.quadrature_point(q)
                      << ": "
                      << "FEFaceValues says " << normal_feval
                      << " but manifold says " << normal_manifold << std::endl;
          }
      }
  deallog << "OK" << std::endl;

  return 0;
}
