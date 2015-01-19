// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2014 by the deal.II authors
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


// make sure that we write boundary lines marked with a non-zero boundary
// indicator correctly in MSH format

#include "../tests.h"
#include <deal.II/base/geometry_info.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <iomanip>


std::ofstream logfile("output");


template <int dim>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  tria.begin_active()->line(0)->set_boundary_indicator(1);
  tria.begin_active()->face(2*dim-1)->set_boundary_indicator(2);

  GridOut grid_out;
  GridOutFlags::Msh flags;
  flags.write_lines = flags.write_faces = true;
  grid_out.set_flags (flags);
  grid_out.write_msh (tria, logfile);
}


int main ()
{
  deallog << std::setprecision (2);
  logfile << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2> ();
  test<3> ();
}

