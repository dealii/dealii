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

// Check that gmsh api correctly reads and writes a mesh with manifold
// information, in all coordinate dimensions

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  deallog << "Testing hypercube in dimensions " << '<' << dim << ',' << spacedim
          << '>' << std::endl;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria, 0, 1, true);
  tria.begin()->face(0)->set_manifold_id(1);

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

  // Generate and print all hypercubes in all dimension combinations
  test<1, 1>();
  test<1, 2>();
  test<1, 3>();

  test<2, 2>();
  test<2, 3>();

  test<3, 3>();
}
