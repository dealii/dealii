// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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

// Compare SphericalManifold::normal_vector with the result of
// FEFaceValues::normal_vector. Those should give very similar results, but an
// old implementation in the manifold gave wrong points

#include <deal.II/base/point.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_generic.h>

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
  FlatManifold<3>        flat;

  Triangulation<dim> tria;
  GridGenerator::half_hyper_shell(tria, Point<dim>(), .5, 1., 0, true);
  tria.set_manifold(0, spherical);
  tria.set_manifold(1, spherical);
  tria.set_manifold(2, spherical);

  tria.refine_global(1);

  MappingQGeneric<dim>   mapping(4);
  QGaussLobatto<dim - 1> quadrature(4);
  FE_Nothing<dim>        dummy;
  FEFaceValues<dim>      fe_values(mapping,
                              dummy,
                              quadrature,
                              update_normal_vectors | update_quadrature_points);

  for (auto &cell : tria.active_cell_iterators())
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->at_boundary(f))
        {
          fe_values.reinit(cell, f);

          // all points should coincide to an accuracy of at least 1e-7 (note
          // that we use a 4-th degree mapping, so its accuracy on the once
          // refined version of the sphere should be enough)
          const double tolerance = 1e-7;
          for (unsigned int q = 0; q < quadrature.size(); ++q)
            {
              const Tensor<1, dim> normal_manifold =
                spherical.normal_vector(cell->face(f),
                                        fe_values.quadrature_point(q));
              if (cell->face(f)->boundary_id() == 2)
                {
                  const Tensor<1, dim> normal_flat =
                    flat.normal_vector(cell->face(f),
                                       fe_values.quadrature_point(q));
                  if (std::abs(1.0 - std::abs(normal_manifold * normal_flat)) >
                      tolerance)
                    deallog
                      << "Error in point " << fe_values.quadrature_point(q)
                      << ": "
                      << "FlatManifold says " << normal_flat
                      << " but SphericalManifold says " << normal_manifold
                      << std::endl
                      << "error: "
                      << std::abs(1.0 - std::abs(normal_manifold * normal_flat))
                      << std::endl;
                }
              else
                {
                  const Tensor<1, dim> normal_feval =
                    fe_values.normal_vector(q);
                  if (std::abs(1.0 - std::abs(normal_manifold * normal_feval)) >
                      tolerance)
                    deallog
                      << "Error in point " << fe_values.quadrature_point(q)
                      << ": "
                      << "FEFaceValues says " << normal_feval
                      << " but manifold says " << normal_manifold << std::endl
                      << "error: "
                      << std::abs(1.0 -
                                  std::abs(normal_manifold * normal_feval))
                      << std::endl;
                }
            }
        }
  deallog << "OK" << std::endl;

  return 0;
}
