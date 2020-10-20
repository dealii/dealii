// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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


// This test verifies that copying a triangulation with an attached
// transfinite interpolation works

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

using namespace dealii;

int
main()
{
  initlog();
  deallog << std::setprecision(9);

  constexpr unsigned int dim = 3;

  dealii::Triangulation<dim> tria_0;
  dealii::Triangulation<dim> tria_1;

  dealii::GridGenerator::torus(tria_0, 2.0, 0.5);
  tria_0.refine_global(1);

  const auto manifold_ids = tria_0.get_manifold_ids();
  for (const auto manifold_id : manifold_ids)
    if (manifold_id != dealii::numbers::flat_manifold_id)
      {
        auto manifold = tria_0.get_manifold(manifold_id).clone();
        tria_1.set_manifold(manifold_id, *manifold);
      }

  const auto construction_data = dealii::TriangulationDescription::Utilities::
    create_description_from_triangulation(tria_0, MPI_COMM_SELF);
  tria_1.create_triangulation(construction_data);

  tria_1.refine_global(1);

  return 0;
}
