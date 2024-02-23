// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>

#include "../tests.h"

void
check_remove_anisotropy()
{
  Point<3>(corners)[3];

  corners[0] = Point<3>(1, 0, 0);
  corners[1] = Point<3>(0, 4, 0);
  corners[2] = Point<3>(0, 0, 2);

  const unsigned int n_subdivisions = 1;

  Triangulation<3> triangulation;
  GridGenerator::subdivided_parallelepiped(triangulation,
                                           n_subdivisions,
                                           corners);
  dealii::GridTools::remove_anisotropy<3>(triangulation,
                                          /*max ratio =*/1.2);

  GridOut grid_out;
  grid_out.write_vtk(triangulation, deallog.get_file_stream());

  triangulation.clear();
}

int
main()
{
  initlog();
  check_remove_anisotropy();
}
