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
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(9);

  constexpr unsigned int dim = 2;

  Triangulation<dim> coarse_triangulation;
  GridGenerator::hyper_cube(coarse_triangulation);

  coarse_triangulation.set_all_manifold_ids(9);
  TransfiniteInterpolationManifold<2> transfinite;
  transfinite.initialize(coarse_triangulation);
  coarse_triangulation.set_manifold(9, transfinite);

  Triangulation<dim> triangulation;
  triangulation.copy_triangulation(coarse_triangulation);

  triangulation.refine_global(2);
}
