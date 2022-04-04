// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2022 by the deal.II authors
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

// Test loading a .msh file from gmsh with some triangular cells having
// negative volume

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <fstream>
#include <iostream>
#include <map>

#include "../tests.h"

using namespace dealii;


template <int dim>
void
print_mesh_info(const Triangulation<dim> &triangulation,
                const std::string &       filename)
{
  deallog << "Mesh info:" << std::endl
          << " dimension: " << dim << std::endl
          << " no. of cells: " << triangulation.n_active_cells() << std::endl;

  {
    std::map<types::boundary_id, unsigned int> boundary_count;
    for (const auto &face : triangulation.active_face_iterators())
      if (face->at_boundary())
        boundary_count[face->boundary_id()]++;

    deallog << " boundary indicators: ";
    for (const std::pair<const types::boundary_id, unsigned int> &pair :
         boundary_count)
      {
        deallog << pair.first << "(" << pair.second << " times) ";
      }
    deallog << std::endl;
  }

  // Finally, produce a graphical representation of the mesh to an output
  // file:
  std::ofstream out(filename);
  GridOut       grid_out;
  grid_out.write_vtu(triangulation, out);
  deallog << " written to " << filename << std::endl << std::endl;
}


void
grid()
{
  Triangulation<2> triangulation;

  GridIn<2> gridin;
  gridin.attach_triangulation(triangulation);
  std::ifstream f(SOURCE_DIR "/grid_in_gmsh_03.msh");
  gridin.read_msh(f);

  print_mesh_info(triangulation, "grid.vtu");
}


int
main()
{
  initlog();
  grid();
}
