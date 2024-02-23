// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Consider a Triangulation with a connected DoFHandler.
// When Triangulation::create_triangulation() gets called (e.g. in
// GridGenerator::merge_triangulations()) the corresponding create signal
// gets triggered.
// At some point, we forgot to reset the DoFHandler at this stage.


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  DoFHandler<dim>    dof_handler(triangulation);

  /*Make a square*/
  Point<dim> point_1, point_2;
  point_1[0] = 0;
  point_1[1] = 0;
  point_2[0] = 1;
  point_2[1] = 1;
  GridGenerator::hyper_rectangle(triangulation, point_1, point_2);

  Triangulation<dim> triangulation_temp;
  point_1[0] = 1;
  point_2[0] = 2;
  GridGenerator::hyper_rectangle(triangulation_temp, point_1, point_2);
  /*glue squares together*/
  GridGenerator::merge_triangulations(triangulation_temp,
                                      triangulation,
                                      triangulation);

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  deallog.push("2d");
  test<2>();
  deallog.pop();
}
