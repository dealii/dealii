// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this test is modeled after after the nedelec_1 test, since that one failed. I
// just wanted to check that the nedelec element that has much of the same
// code also works. turns out, all was fine

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);

  FE_Nedelec<dim> fe(0);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  QGauss<dim - 1>   q(2);
  FEFaceValues<dim> fe_values(fe, q, update_values | update_gradients);
  fe_values.reinit(dof_handler.begin_active(), 0);

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test<2>();
  test<3>();

  return 0;
}
