// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2011 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

// Check if Triangulation<1>::get_boundary_ids() works for 1d grids.


int
main()
{
  initlog();

  Triangulation<1> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);
  const std::vector<types::boundary_id> indicators =
    triangulation.get_boundary_ids();
  for (unsigned int i = 0; i < indicators.size(); ++i)
    deallog << int(indicators[i]) << std::endl;

  return 0;
}
