// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check SphericalManifold::normal_vector on a set of points of a
// hypersphere mesh

#include <deal.II/base/point.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  constexpr unsigned int dim = 3;
  SphericalManifold<3>   spherical;

  Triangulation<dim> tria;
  GridGenerator::hyper_shell(tria, Point<dim>(), 0.5, 1., 96, true);
  tria.set_all_manifold_ids(0);
  tria.set_manifold(0, spherical);

  for (auto &cell : tria.active_cell_iterators())
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->at_boundary(f))
        {
          if (std::abs(cell->face(f)->vertex(1).norm() - 1.) < 1e-1)
            {
              // at outer boundary, the normal should point outwards and be
              // aligned with the point down to roundoff. Check this for the
              // second face vertex as well as the center.
              deallog
                << "Outer normal correctness (should be 0): "
                << (1.0 - cell->face(f)->vertex(1) *
                            spherical.normal_vector(cell->face(f),
                                                    cell->face(f)->vertex(1)))
                << ' '
                << (1.0 -
                    cell->face(f)->center(/*respect_manifold=*/true) *
                      spherical.normal_vector(cell->face(f),
                                              cell->face(f)->center(true)))
                << std::endl;
            }
          else if (std::abs(cell->face(f)->vertex(1).norm() - 0.5) < 1e-1)
            {
              // at inner boundary, the normal should point inwards and be
              // aligned with the point down to roundoff. Check this for the
              // second face vertex as well as the center.
              deallog
                << "Inner normal correctness (should be 0): "
                << (0.5 - cell->face(f)->vertex(1) *
                            spherical.normal_vector(cell->face(f),
                                                    cell->face(f)->vertex(1)))
                << ' '
                << (0.5 -
                    cell->face(f)->center(/*respect_manifold=*/true) *
                      spherical.normal_vector(cell->face(f),
                                              cell->face(f)->center(true)))
                << std::endl;
            }
        }
  for (auto &cell : tria.active_cell_iterators())
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->at_boundary(f))
        {
          // the first choice will contain perturbations in the normal because
          // the center without respecting the manifold is not on the face to
          // which the algorithm computes its points. As an alternative, we
          // adjust the point to lie on the same radius as the vertex points
          // by a weighting where the normal should again be perfectly aligned
          // with the point itself.
          deallog
            << "Approximate normal in " << cell->face(f)->center(false) << ":  "
            << spherical.normal_vector(cell->face(f),
                                       cell->face(f)->center(false))
            << ", adjusted: "
            << spherical.normal_vector(cell->face(f),
                                       cell->face(f)->center(false) /
                                         cell->face(f)->center(false).norm() *
                                         cell->face(f)->vertex(1).norm())
            << std::endl;
        }

  return 0;
}
