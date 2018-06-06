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
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <fstream>

using namespace dealii;

int
main()
{
  Triangulation<2, 3> tria;
  GridIn<2, 3>        gridin;
  gridin.attach_triangulation(tria);

  {
    std::ofstream ofile("plane.obj");
    ofile << "mtllib plane.mtl" << std::endl
          << "o Plane" << std::endl
          << "v -1.000000 0.000000 1.000000" << std::endl
          << "v 1.000000 0.000000 1.000000" << std::endl
          << "v -1.000000 0.000000 -1.000000" << std::endl
          << "v 1.000000 0.000000 -1.000000" << std::endl
          << "vn 0.0000 1.0000 0.0000" << std::endl
          << "usemtl None" << std::endl
          << "s off" << std::endl
          << "f 1//1 2//1 4//1 3//1" << std::endl;
  }
  gridin.read_assimp("plane.obj", -1, true, 1e-3);
  return 0;
}
