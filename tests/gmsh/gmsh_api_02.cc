// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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

// Create all reference cells, and output them with gmsh.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test(const std::uint8_t kind, const std::string out = "")
{
  deallog << "Testing kind(" << (int)kind << ") in dimensions " << '<' << dim
          << ',' << spacedim << '>' << std::endl;

  Triangulation<dim, spacedim> tria;
  GridGenerator::reference_cell(
    tria, internal::ReferenceCell::make_reference_cell_from_int(kind));
  GridOut go;
  go.write_msh(tria, "output.msh");
  cat_file("output.msh");
  // Uncomment to inspect each file separately
  // if (out != "")
  //   go.write_msh(tria, out);
  (void)out;
}

int
main()
{
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
