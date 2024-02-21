// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// The first attempt to fix the _01 test had a bug in it in that we
// managed to invalidate the cell stored in FEValues but forgot to
// disconnect from the triangulation's signal. We were therefore still
// getting reminders from the triangulation even though we weren't
// doing anything with it any more

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  FE_Q<dim>         fe(1);
  const QGauss<dim> quadrature(2);
  FEValues<dim>     fe_values(fe, quadrature, update_values);

  // initialize FEValues with the first cell
  fe_values.reinit(tr.begin_active());

  // then invalidate the cell iterator
  tr.refine_global(1);
  // and invalidate it again. this
  // shouldn't do any further harm
  tr.refine_global(1);

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test<1>();
  test<2>();
  test<3>();
}
