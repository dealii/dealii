// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2022 by the deal.II authors
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


// This test verifies that the transfinite interpolation can be initialized
// in the codimension one case.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <memory>

#include "../tests.h"

int
main()
{
  initlog();

  using namespace dealii;

  Triangulation<2, 3> tria;

  GridGenerator::hyper_cube(tria, -1, 1);
  TransfiniteInterpolationManifold<2, 3> transfinite;
  transfinite.initialize(tria);

  deallog << "OK" << std::endl;

  return 0;
}
