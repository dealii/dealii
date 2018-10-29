// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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
