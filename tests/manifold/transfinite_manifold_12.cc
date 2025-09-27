// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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


  Triangulation<2, 3> tria;

  GridGenerator::hyper_cube(tria, -1, 1);
  TransfiniteInterpolationManifold<2, 3> transfinite;
  transfinite.initialize(tria);

  deallog << "OK" << std::endl;

  return 0;
}
