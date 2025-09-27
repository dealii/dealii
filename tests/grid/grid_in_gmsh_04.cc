// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test that we can read dim = 1, codim = 2 meshes with boundary ids correctly.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <fstream>
#include <iostream>
#include <map>

#include "../tests.h"

template <int dim, int spacedim = dim>
void
print_mesh_info(const Triangulation<dim, spacedim> &triangulation,
                const std::string                  &filename)
{
  deallog << "Mesh info:" << std::endl
          << " dimension: " << dim << std::endl
          << " no. of cells: " << triangulation.n_active_cells() << std::endl;

  {
    std::map<types::boundary_id, unsigned int> boundary_count;
    for (const auto &cell : triangulation.active_cell_iterators())
      for (const unsigned int face_no : cell->face_indices())
        if (cell->face(face_no)->at_boundary())
          boundary_count[cell->face(face_no)->boundary_id()]++;

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


int
main()
{
  initlog();
  std::ifstream in(SOURCE_DIR "/grid_in_gmsh_04.msh");

  Triangulation<1, 2> tria;
  GridIn<1, 2>        gi;
  gi.attach_triangulation(tria);
  gi.read_msh(in);

  print_mesh_info(tria, "grid_in_gmsh_04.vtu");
}
