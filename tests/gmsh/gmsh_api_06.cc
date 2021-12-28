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

// Create hyper cube, 0bcs, write, read, and write again with gmsh

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
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
  GridGenerator::hyper_cube(tria, 0, 1, true);

  GridOut go;
  go.write_msh(tria, "output.msh");

  Triangulation<dim, spacedim> tria2;
  GridIn<dim, spacedim>        gi(tria2);
  gi.read_msh("output.msh");

  go.write_msh(tria2, "output2.msh");

  deallog << "Original mesh: " << std::endl;
  cat_file("output.msh");
  deallog << "Regenerated mesh: " << std::endl;
  cat_file("output2.msh");
}