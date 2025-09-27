// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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
