// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// take a 3d mesh, take all vertices, shift them a little bit and check that
// we correctly identify the closest vertex position
// The result should be an increasing sequence of numbers

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



void
check(Triangulation<3> &tria)
{
  const std::vector<Point<3>> &v = tria.get_vertices();
  for (unsigned i = 0; i < v.size(); ++i)
    deallog << '['
            << GridTools::find_closest_vertex(tria,
                                              v[i] +
                                                Point<3>(0.01, -0.01, 0.01))
            << "] ";

  deallog << std::endl;
}


int
main()
{
  initlog();

  try
    {
      Triangulation<3> coarse_grid;
      GridGenerator::hyper_cube(coarse_grid);
      coarse_grid.refine_global(3);
      check(coarse_grid);
    }
  catch (const std::exception &exc)
    {
      // we shouldn't get here...
      deallog << "Caught an error..." << std::endl;
      deallog << exc.what() << std::endl;
    }
}
