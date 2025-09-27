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

// Create a hyper ball, use gmsh to output it and to read it back in.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test(const Triangulation<dim, spacedim> &tria)
{
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



int
main(int argc, char **argv)
{
  // gmsh might be build with mpi support enabled.
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  initlog();
  {
    Triangulation<2> tria;
    GridGenerator::hyper_ball(tria);
    test(tria);
  }
  {
    Triangulation<2, 3> tria;
    GridGenerator::hyper_sphere(tria);
    test(tria);
  }
  {
    Triangulation<3, 3> tria;
    GridGenerator::hyper_ball(tria);
    test(tria);
  }
}
