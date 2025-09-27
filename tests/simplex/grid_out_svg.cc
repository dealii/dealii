// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Write a file in the SVG format with (linear) triangular and tetrahedral
// elements created by the GMSH program.

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <fstream>

#include "../tests.h"

template <int dim, int spacedim>
void
check_file(const std::string &file_name)
{
  Triangulation<dim, spacedim> tria;

  GridIn<dim, spacedim> grid_in;
  grid_in.attach_triangulation(tria);
  std::ifstream input_file(file_name);
  grid_in.read_vtk(input_file);

  GridOut grid_out;
#if false
  std::ofstream out("mesh-" + std::to_string(spacedim) + ".svg");
  grid_out.write_svg(tria, out);
#else
  grid_out.write_svg(tria, deallog.get_file_stream());
#endif

  deallog << "OK!" << std::endl;
}


int
main()
{
  initlog();
  // TRIANGULAR ELEMENTS
  // dim = spacedim = 2
  deallog.push("triangluar_elements_dim2_spacedim2: ");
  check_file<2, 2>(std::string(SOURCE_DIR "/grid_in_vtk/tri.vtk"));
  deallog.pop();
}
