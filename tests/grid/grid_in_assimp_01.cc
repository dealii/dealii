//-----------------------------------------------------------
//
//    Copyright (C) 2017 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE at
//    the top level of the deal.II distribution.
//
//-----------------------------------------------------------

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


int
main()
{
  initlog();
  Triangulation<2, 3> tria;
  GridIn<2, 3>        gridin;
  gridin.attach_triangulation(tria);
  gridin.read_assimp(SOURCE_DIR "/grids/torus.obj", -1, true, 1e-3);

  GridOut go;
  go.write_msh(tria, deallog.get_file_stream());
  return 0;
}
