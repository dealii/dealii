// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
