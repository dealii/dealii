// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test that transfinite interpolation manifold returns valid results on the
// same geometry as transfinite_manifold_05 (ball inside square, curved
// spherical surface). In an initial version, the line search in Newton for
// the transfinite interpolation would eagerly search too far outside the
// valid chart domain, leading to failures in the spherical manifold.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"


int
main()
{
  initlog();

  const int          dim = 3;
  Triangulation<dim> tria1, tria2, tria;
  GridGenerator::hyper_shell(tria1, Point<dim>(), 0.4, std::sqrt(dim), 6);
  GridGenerator::hyper_ball(tria2, Point<dim>(), 0.4);
  GridGenerator::merge_triangulations(tria1, tria2, tria);
  tria.set_all_manifold_ids(0);
  for (typename Triangulation<dim>::cell_iterator cell = tria.begin();
       cell != tria.end();
       ++cell)
    {
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        {
          bool face_at_sphere_boundary = true;
          for (unsigned int v = 0; v < GeometryInfo<dim - 1>::vertices_per_cell;
               ++v)
            if (std::abs(cell->face(f)->vertex(v).norm() - 0.4) > 1e-12)
              face_at_sphere_boundary = false;
          if (face_at_sphere_boundary)
            cell->face(f)->set_all_manifold_ids(1);
        }
    }
  static const SphericalManifold<dim> spherical_manifold;
  tria.set_manifold(1, spherical_manifold);
  static TransfiniteInterpolationManifold<dim> transfinite;
  transfinite.initialize(tria);
  tria.set_manifold(0, transfinite);

  tria.refine_global(1);

  deallog.precision(10);
  deallog << "Cell centers" << std::endl;
  for (auto &cell : tria.cell_iterators())
    deallog << cell->id() << " has center "
            << cell->center(/*respect_manifold*/ true) << std::endl;

  return 0;
}
