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

// Create all reference cells, and output them with gmsh.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test(const std::uint8_t kind, const std::string out = "")
{
  deallog << "Testing kind(" << (int)kind << ") in dimensions " << '<' << dim
          << ',' << spacedim << '>' << std::endl;

  Triangulation<dim, spacedim> tria;
  GridGenerator::reference_cell(tria,
                                internal::make_reference_cell_from_int(kind));

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

  // Uncomment to inspect each file separately
  // if (out != "")
  //   go.write_msh(tria2, out);
  (void)out;
}

int
main(int argc, char **argv)
{
  // gmsh might be build with mpi support enabled.
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  initlog();

  // Generate and print all reference cells
  // Lines
  test<1, 1>(1, "line_11.msh");
  test<1, 2>(1, "line_12.msh");
  test<1, 3>(1, "line_13.msh");

  // Triangles
  test<2, 2>(2, "tri_22.msh");
  test<2, 3>(2, "tri_23.msh");

  // Quads
  test<2, 2>(3, "quad_22.msh");
  test<2, 3>(3, "quad_23.msh");

  // Tet
  test<3, 3>(4, "tet_3.msh");

  // Pyramid
  test<3, 3>(5, "pyramid_3.msh");

  // Wedge
  test<3, 3>(6, "wedge_3.msh");

  // Hex
  test<3, 3>(7, "hex_3.msh");
}
