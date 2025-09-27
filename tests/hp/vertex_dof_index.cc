// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check cell->vertex_dof_index() in the hp case when called on a cell
// where the active_fe_index should be clear.


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test()
{
  // create 2 triangulations with the
  // same coarse grid, and refine
  // them differently
  Triangulation<dim> tria;

  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  hp::FECollection<dim> fe;
  fe.push_back(FE_Q<dim>(1));
  fe.push_back(FE_Q<dim>(2));
  DoFHandler<dim> dh(tria);
  dh.begin_active()->set_active_fe_index(1);
  dh.distribute_dofs(fe);

  deallog << dh.begin_active()->vertex_dof_index(0, 0) << std::endl;
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
