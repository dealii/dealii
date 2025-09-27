// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test()
{
  // create 3 triangulations
  Triangulation<dim> tria[3];

  GridGenerator::hyper_cube(tria[0]);
  tria[0].refine_global(1);

  GridGenerator::hyper_cube(tria[1]);
  GridTools::scale(2, tria[1]);
  tria[1].refine_global(2);

  if (dim != 1)
    GridGenerator::hyper_ball(tria[2]);
  else
    {
      GridGenerator::hyper_cube(tria[2]);
      GridTools::shift(Point<dim>(2.), tria[2]);
    }

  tria[2].refine_global(3);

  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      {
        Assert(GridTools::have_same_coarse_mesh(tria[i], tria[j]) == (i == j),
               ExcInternalError());

        deallog << "meshes " << i << " and " << j << ": "
                << (GridTools::have_same_coarse_mesh(tria[i], tria[j]) ?
                      "true" :
                      "false")
                << std::endl;
      }
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
