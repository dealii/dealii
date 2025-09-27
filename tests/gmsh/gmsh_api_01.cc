// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Create a hyper ball, and use gmsh to output it.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

int
main(int argc, char **argv)
{
  // gmsh might be build with mpi support enabled.
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
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
