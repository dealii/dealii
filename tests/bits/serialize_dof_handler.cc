// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2018 by the deal.II authors
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



// Test an issue reported by David Wells: serializing a DoFHandler object did
// not work right out of the box without manually including additional header
// files

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <boost/archive/text_oarchive.hpp>

#include <iostream>
#include <string>

#include "../tests.h"


int
main()
{
  initlog();

  Triangulation<2> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(4);
  DoFHandler<2> dof_handler(triangulation);
  FE_Q<2>       finite_element(1);
  dof_handler.distribute_dofs(finite_element);

  std::ostringstream            out_stream;
  boost::archive::text_oarchive archive(out_stream);

  archive << dof_handler;
  dof_handler.clear();

  deallog << "OK" << std::endl;
}
