// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check PolarManifold for normal vector issues. The PolarManifold
// should produce the same normal vector as the SphericalManifold in 2D,
// in particular for points at the same radius.

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"


int
main()
{
  initlog();

  // Center and radius of the Ball

  const PolarManifold<2, 2>     polar_manifold;
  const SphericalManifold<2, 2> spherical_manifold;

  Triangulation<2, 2> tria;
  GridGenerator::hyper_shell(tria, Point<2>(), .3, .6, 12);

  tria.set_manifold(1, polar_manifold);

  for (typename Triangulation<2, 2>::active_cell_iterator cell =
         tria.begin_active();
       cell != tria.end();
       ++cell)
    {
      cell->set_all_manifold_ids(1);
    }

  for (typename Triangulation<2, 2>::active_cell_iterator cell =
         tria.begin_active();
       cell != tria.end();
       ++cell)
    {
      for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
        {
          const Point<2>     position = cell->face(f)->center();
          const Tensor<1, 2> polar_normal =
            polar_manifold.normal_vector(cell->face(f), position);
          const Tensor<1, 2> spherical_normal =
            spherical_manifold.normal_vector(cell->face(f), position);
          deallog << "Position: " << position << std::endl;
          deallog << "Polar normal: " << polar_normal << std::endl;
          deallog << "Spherical normal: " << spherical_normal << std::endl;
          deallog << "Relative difference: "
                  << (polar_normal - spherical_normal).norm() /
                       polar_normal.norm()
                  << std::endl
                  << std::endl;
        }
    }

  return 0;
}
