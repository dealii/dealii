// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
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



// Test an issue reported by David Wells: serializing a DoFHandler object did
// not work right out of the box without manually including additional header
// files

#include "../tests.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <boost/archive/text_oarchive.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <fstream>


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Triangulation<2> triangulation;
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (4);
  DoFHandler<2> dof_handler (triangulation);
  FE_Q<2> finite_element (1);
  dof_handler.distribute_dofs (finite_element);

  std::ostringstream out_stream;
  boost::archive::text_oarchive archive(out_stream);

  archive << dof_handler;
  dof_handler.clear ();

  deallog << "OK" << std::endl;
}
