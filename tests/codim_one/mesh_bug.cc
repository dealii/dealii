// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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
    if (vertices[i](2) > 1e-7)
      std::cout << "Error!" << std::endl;


  GridOut go;
  go.write_ucd(tria, deallog.get_file_stream());

  return 0;
}
