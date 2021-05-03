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

// Create a hyper ball, and use gmsh to output it.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

int
main()
{
  initlog();

  const unsigned int dim      = 2;
  const unsigned int spacedim = 2;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_ball(tria);
  tria.refine_global(1);

  GridOut go;
  go.write_msh(tria, "output.msh");

  cat_file("output.msh");
}
