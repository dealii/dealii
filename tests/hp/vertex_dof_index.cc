// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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
