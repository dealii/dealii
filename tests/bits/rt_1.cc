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



// there was a bug in the RT element that Oliver Kayser-Herold fixed in
// January 2005. Check this

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


template <int dim>
void
test(const unsigned int degree, const unsigned int q_order)
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);

  FE_RaviartThomas<dim> fe(degree);
  DoFHandler<dim>       dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  QGauss<dim - 1>   q(q_order);
  FEFaceValues<dim> fe_values(fe,
                              q,
                              update_values | update_gradients |
                                update_hessians | update_quadrature_points |
                                update_jacobians);
  fe_values.reinit(dof_handler.begin_active(), 0);

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  for (unsigned int degree = 0; degree < 3; ++degree)
    for (unsigned int q_order = 1; q_order <= 3; ++q_order)
      test<2>(degree, q_order);

  return 0;
}
