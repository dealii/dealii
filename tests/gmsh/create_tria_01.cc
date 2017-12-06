//-----------------------------------------------------------
//
//    Copyright (C) 2014 - 2017 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

// Read a file in iges format, and write it out again in the same
// format.

#include "../tests.h"

int main ()
{
  initlog();
  std::ofstream geo("file.geo");

  geo << "Lx = 25.0;" << std::endl
      << "Ly = 1.0;" << std::endl
      << "Point(1) = {0, 0, 0, 1};" << std::endl
      << "Point(2) = {Lx, 0, 0, 1};" << std::endl
      << "Point(3) = {Lx, Ly, 0, 1};" << std::endl
      << "Point(4) = {0, Ly, 0, 1};" << std::endl
      << "Line(1) = {1, 2};" << std::endl
      << "Line(2) = {2, 3};" << std::endl
      << "Line(3) = {3, 4};" << std::endl
      << "Line(4) = {4, 1};" << std::endl
      << "Line Loop(5) = {3, 4, 1, 2};" << std::endl
      << "Plane Surface(0) = {5};" << std::endl
      << "Transfinite Surface{0} = {5};" << std::endl
      << "Recombine Surface {0};" << std::endl;

  geo.close();

  std::system(DEAL_II_GMSH_EXECUTABLE_PATH " -2 file.geo 1>file.log 2>file_warn.log");

  cat_file("file.msh");

  return 0;
}
