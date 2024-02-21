// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// find_active_cell_around_point should throw an exception if the
// point is outside. Test that.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <time.h>

#include <iostream>
#include <list>
#include <sstream>
#include <string>

#include "../tests.h"


void
test()
{
  Triangulation<2> tr;
  GridGenerator::hyper_cube(tr);

  Point<2> p;
  p[0] = -0.1;
  p[1] = 0.5;

  MappingQ<2> mapping(1);

  auto c = GridTools::find_active_cell_around_point(mapping, tr, p);
  if (c.first.state() != IteratorState::valid)
    deallog << "outside" << std::endl;
  deallog << "done" << std::endl;
}

int
main(int argc, char **argv)
{
  initlog();

  test();

  return 0;
}
