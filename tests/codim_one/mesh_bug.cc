// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// a short (a few lines) description of what the program does

#include "../tests.h"

// all include files you need here


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>


int
main()
{
  initlog();

  Triangulation<2, 3> tria;
  GridIn<2, 3>        gi;
  gi.attach_triangulation(tria);
  std::ifstream infile(SOURCE_DIR "/mesh_bug.output");
  gi.read(infile);

  const std::vector<Point<3>> &vertices = tria.get_vertices();

  for (unsigned int i = 0; i < vertices.size(); ++i)
    if (vertices[i][2] > 1e-7)
      std::cout << "Error!" << std::endl;


  GridOut go;
  go.write_ucd(tria, deallog.get_file_stream());

  return 0;
}
