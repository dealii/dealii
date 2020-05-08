// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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
