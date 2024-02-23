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


// Read in a dim=1, spacedim=3 file in VTK format, and replicate
// it. This failed at one point because we erroneously tried to find
// correct directions for each line.


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

  Triangulation<1, 3> triangulation2;
  GridIn<1, 3>        grid_in1;
  grid_in1.attach_triangulation(triangulation2);
  std::ifstream input_file1(SOURCE_DIR "/grid_in_vtk_1d_3d.tria2.vtk");
  grid_in1.read_vtk(input_file1);

  Triangulation<1, 3> triangulation3;
  GridGenerator::replicate_triangulation(triangulation2, {2}, triangulation3);

  GridOut go;
  go.write_gnuplot(triangulation3, deallog.get_file_stream());
}
