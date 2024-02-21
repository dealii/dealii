// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test that the GMSH executable actually works.
// This is a reduced version of gmsh/create_tria_01.

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>

#include <fstream>

int
main()
{
  std::ofstream geo("file.geo");

  geo << "Lx = 25.0;" << std::endl
      << "Ly = 1.0;" << std::endl
      << "Point(1) = {0, 0, 0, Lx};" << std::endl
      << "Point(2) = {Lx, 0, 0, Lx};" << std::endl
      << "Point(3) = {Lx, Ly, 0, Lx};" << std::endl
      << "Point(4) = {0, Ly, 0, Lx};" << std::endl
      << "Line(1) = {1, 2};" << std::endl
      << "Line(2) = {2, 3};" << std::endl
      << "Line(3) = {3, 4};" << std::endl
      << "Line(4) = {4, 1};" << std::endl
      << "Transfinite Line {1} = 2 Using Progression 1;" << std::endl
      << "Transfinite Line {2} = 2 Using Progression 1;" << std::endl
      << "Transfinite Line {3} = 2 Using Progression 1;" << std::endl
      << "Transfinite Line {4} = 2 Using Progression 1;" << std::endl
      << "Line Loop(5) = {3, 4, 1, 2};" << std::endl
      << "Plane Surface(6) = {5};" << std::endl
      << "Transfinite Surface {6};" << std::endl
      << "Recombine Surface {6};" << std::endl
      << "Physical Surface(\"dealii_surface\") = {6};" << std::endl;

  geo.close();

  const int ierr = std::system(DEAL_II_GMSH_EXECUTABLE_PATH
                               " -2 file.geo 1>file.log 2>file_warn.log");
  (void)ierr;
  Assert(ierr == 0, dealii::ExcInternalError());

  std::remove("file.geo");
  std::remove("file.log");
  std::remove("file_warn.log");
  return 0;
}
