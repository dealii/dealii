// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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


#include "../tests.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <iomanip>
#include <string>

std::ofstream logfile("output");


template <int dim>
void test (const std::string &infilename)
{
  Triangulation<dim> tria;
  GridIn<dim> gi;
  gi.attach_triangulation (tria);
  gi.read (infilename);

  logfile<<"------------------------------------------"<<std::endl;

  GridOut grid_out;
  grid_out.set_flags (GridOutFlags::Ucd(true));
  grid_out.write_ucd (tria, logfile);
}

int main ()
{
  test<2> (std::string(SOURCE_DIR "/grid_in_tecplot/1.dat"));
  test<2> (std::string(SOURCE_DIR "/grid_in_tecplot/2.dat"));
  test<2> (std::string(SOURCE_DIR "/grid_in_tecplot/3.dat"));
  test<2> (std::string(SOURCE_DIR "/grid_in_tecplot/4.dat"));
}

