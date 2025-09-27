// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2025 by the deal.II authors
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
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



void
test()
{
  const int dim      = 2;
  const int spacedim = 3;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_ball(tria, Point<spacedim>(1, 1, 1), 0.5);
  tria.refine_global(1);

  GridOut grid_out;
  grid_out.write_gnuplot(tria, deallog.get_file_stream());
}


int
main()
{
  initlog();
  test();
}
