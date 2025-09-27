// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test loading a .msh file from gmsh with some cells having negative volume

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



template <int dim>
void
print_mesh_info(const Triangulation<dim> &triangulation,
                const std::string        &filename)
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
  std::ifstream f(SOURCE_DIR "/grid_in_gmsh_02.msh");
  gridin.read_msh(f);

  print_mesh_info(triangulation, "grid.vtu");
}


int
main()
{
  initlog();
  grid();
}
