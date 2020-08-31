// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// Write a file in the EPS format with (linear) triangular and tetrahedral
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
  std::ofstream out("mesh-" + std::to_string(spacedim) + ".eps");
  grid_out.write_eps(tria, out);
#else
  grid_out.write_eps(tria, deallog.get_file_stream());
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

  // TETRAHEDRAL ELEMENTS
  // dim = spacedim = 3
  deallog.push("tetrahedral_elements_dim3_spacedim3: ");
  check_file<3, 3>(std::string(SOURCE_DIR "/grid_in_vtk/tet.vtk"));
  deallog.pop();
}
