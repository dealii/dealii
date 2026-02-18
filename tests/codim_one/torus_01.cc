// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2008 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// create a torus mesh and refine it.

#include "../tests.h"

// all include files you need here

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <string>

int
main()
{
  const int dim      = 2;
  const int spacedim = 3;

  initlog();

  Triangulation<dim, spacedim> tria;
  GridGenerator::torus(tria, 1.5, .5);
  tria.refine_global(2);

  GridOut grid_out;
  grid_out.write_gnuplot(tria, deallog.get_file_stream());

  return 0;
}
