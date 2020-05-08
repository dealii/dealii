// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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


// Read a file in iges format, and write it out again in the same
// format.

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
gmsh_grid(const char *name)
{
  Triangulation<dim> tria;
  GridIn<dim>        grid_in;
  grid_in.attach_triangulation(tria);
  std::ifstream input_file(name);
  grid_in.read_msh(input_file);

  deallog << "  " << tria.n_active_cells() << " active cells" << std::endl;

  int hash  = 0;
  int index = 0;
  for (typename Triangulation<dim>::active_cell_iterator c =
         tria.begin_active();
       c != tria.end();
       ++c, ++index)
    for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
      hash += (index * i * c->vertex_index(i)) % (tria.n_active_cells() + 1);
  deallog << "  hash=" << hash << std::endl;
}

int
main()
{
  initlog();
  std::ofstream geo("file.geo");

  geo << "Lx = 25.0;" << std::endl
      << "Ly = 1.0;" << std::endl
      << "Point(1) = {0, 0, 0, Lx};" << std::endl
      << "Point(2) = {Lx, 0, 0, Lx};" << std::endl
      << "Point(3) = {Lx, Ly, 0, Lx};" << std::endl
      << "Point(4) = {0, Ly, 0, Lx};" << std::endl
      << "Line(1) = {1, 2};" << std::endl
      << "Line(2) = {2, 3};" << std::endl
      << "Line(3) = {3, 4};" << std::endl
      << "Line(4) = {4, 1};" << std::endl
      << "Transfinite Line {1} = 2 Using Progression 1;" << std::endl
      << "Transfinite Line {2} = 2 Using Progression 1;" << std::endl
      << "Transfinite Line {3} = 2 Using Progression 1;" << std::endl
      << "Transfinite Line {4} = 2 Using Progression 1;" << std::endl
      << "Line Loop(5) = {3, 4, 1, 2};" << std::endl
      << "Plane Surface(6) = {5};" << std::endl
      << "Transfinite Surface {6};" << std::endl
      << "Recombine Surface {6};" << std::endl
      << "Physical Surface(\"dealii_surface\") = {6};" << std::endl;

  geo.close();

  const int ierr = std::system(DEAL_II_GMSH_EXECUTABLE_PATH
                               " -2 file.geo 1>file.log 2>file_warn.log");
  Assert(ierr == 0, ExcInternalError());
  gmsh_grid<2>("file.msh");
  cat_file("file.msh");
  std::remove("file.log");
  std::remove("file_warn.log");

  return 0;
}
