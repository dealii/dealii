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

// Test vertex manifold ids in 1D.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

int
main()
{
  initlog(true);

  SphericalManifold<1> spherical_manifold;


  constexpr types::manifold_id spherical_manifold_id = 42;

  Triangulation<1> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);
  triangulation.set_manifold(spherical_manifold_id, spherical_manifold);

  triangulation.begin_active()->set_all_manifold_ids(spherical_manifold_id);

  for (const auto &cell : triangulation.active_cell_iterators())
    {
      deallog << "current cell manifold id: " << cell->manifold_id()
              << std::endl;

      for (const unsigned int vertex_n : GeometryInfo<1>::vertex_indices())
        {
          deallog << "current vertex: " << cell->vertex(vertex_n) << std::endl;
          deallog << "current vertex manifold id: "
                  << cell->face(vertex_n)->manifold_id() << std::endl;
        }
    }
}
