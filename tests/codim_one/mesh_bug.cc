// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// a short (a few lines) description of what the program does

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>

// all include files you need here


#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/grid_out.h>
#include <grid/grid_in.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>


int main ()
{
  std::ofstream logfile("output");

  Triangulation<2,3> tria;
  GridIn<2,3> gi;
  gi.attach_triangulation(tria);
  std::ifstream infile(SOURCE_DIR "/mesh_bug.output");
  gi.read(infile);

  const std::vector<Point<3> > &vertices = tria.get_vertices();

  for (unsigned int i=0; i<vertices.size(); ++i)
    if (vertices[i](2)>1e-7)
      std::cout << "Error!" << std::endl;


  GridOut go;
  go.write_ucd(tria, logfile);

  return 0;
}

