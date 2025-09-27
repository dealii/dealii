// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test1()
{
  // test 1: hypercube
  if (true)
    {
      Triangulation<dim> tria;
      GridGenerator::hyper_cube(tria);

      for (unsigned int i = 0; i < 2; ++i)
        {
          tria.refine_global(2);
          deallog << dim << "d, "
                  << "max diameter: " << GridTools::maximal_cell_diameter(tria)
                  << std::endl;
          Assert(GridTools::maximal_cell_diameter(tria) >=
                   GridTools::minimal_cell_diameter(tria),
                 ExcInternalError());
        };
    };

  // test 2: hyperball
  if (dim >= 2)
    {
      Triangulation<dim> tria;
      GridGenerator::hyper_ball(tria, Point<dim>(), 1);
      tria.reset_manifold(0);

      for (unsigned int i = 0; i < 2; ++i)
        {
          tria.refine_global(2);
          deallog << dim << "d, "
                  << "max diameter: " << GridTools::maximal_cell_diameter(tria)
                  << std::endl;
          Assert(GridTools::maximal_cell_diameter(tria) >=
                   GridTools::minimal_cell_diameter(tria),
                 ExcInternalError());
        };
    };
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test1<1>();
  test1<2>();
  test1<3>();

  return 0;
}
