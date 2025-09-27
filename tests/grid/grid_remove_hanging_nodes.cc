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

#include "../tests.h"

void
check_remove_hanging_nodes()
{
  Point<2> corners[2];

  corners[0]                        = Point<2>(1, 0);
  corners[1]                        = Point<2>(0, 4);
  const unsigned int n_subdivisions = 1;


  Triangulation<2> tria;
  GridGenerator::subdivided_parallelepiped(tria, n_subdivisions, corners);

  tria.refine_global();

  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  tria.refine_global();

  dealii::GridTools::remove_hanging_nodes(tria, /*isotropic=*/false);

  GridOut grid_out;
  grid_out.write_vtk(tria, deallog.get_file_stream());

  tria.clear();
}

int
main()
{
  initlog();
  check_remove_hanging_nodes();
}
