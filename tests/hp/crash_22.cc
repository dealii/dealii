// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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
  point_1(0) = 0;
  point_1(1) = 0;
  point_2(0) = 1;
  point_2(1) = 1;
  GridGenerator::hyper_rectangle(triangulation, point_1, point_2);

  Triangulation<dim> triangulation_temp;
  point_1(0) = 1;
  point_2(0) = 2;
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
