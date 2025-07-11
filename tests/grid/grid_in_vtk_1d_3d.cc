// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Read in two dim=1, spacedim=3 files in VTK format, and merge them
// into one triangulation. This failed at one point because we
// erroneously tried to find correct directions for each line.


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "../tests.h"


int
main()
{
  initlog();

  Triangulation<1, 3> triangulation1;
  GridIn<1, 3>        grid_in1;
  grid_in1.attach_triangulation(triangulation1);
  std::ifstream input_file1(SOURCE_DIR "/grid_in_vtk_1d_3d.tria1.vtk");
  grid_in1.read_vtk(input_file1);

  Triangulation<1, 3> triangulation2;
  GridIn<1, 3>        grid_in2;
  grid_in2.attach_triangulation(triangulation2);
  std::ifstream input_file2(SOURCE_DIR "/grid_in_vtk_1d_3d.tria2.vtk");
  grid_in2.read_vtk(input_file2);

  deallog << "Triangulation 1:\n";
  for (const auto &cell : triangulation1.active_cell_iterators())
    {
      deallog << "  cell=" << cell << std::endl;
      for (const auto c : cell->vertex_indices())
        {
          Point<3> &v = cell->vertex(c);
          deallog << "    vertex: " << v << std::endl;
        }
    }

  for (const auto &boundary_id : triangulation1.get_boundary_ids())
    deallog << "  boundary_id = " << boundary_id << std::endl;

  deallog << "Triangulation 2:\n";
  for (const auto &cell : triangulation2.active_cell_iterators())
    {
      deallog << "  cell=" << cell << std::endl;
      for (const auto c : cell->vertex_indices())
        {
          Point<3> &v = cell->vertex(c);
          deallog << "    vertex: " << v << std::endl;
        }
    }

  for (const auto &boundary_id : triangulation1.get_boundary_ids())
    deallog << "  boundary_id = " << boundary_id << std::endl;

  GridGenerator::merge_triangulations(triangulation2,
                                      triangulation1,
                                      triangulation2);
  GridOut go;
  go.write_gnuplot(triangulation2, deallog.get_file_stream());
}
