// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check that FEValues can be reinit'ed with a FE_RaviartThomasNodal on the
// mesh of a 2d ball when no finite element entries are requested

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include <string>
#include <vector>

#include "../tests.h"

int
main()
{
  initlog();

  const int          dim = 2;
  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria, Point<dim>(), 1.);

  FE_RaviartThomasNodal<dim> fe(0);
  DoFHandler<dim>            dof(tria);
  dof.distribute_dofs(fe);

  FEValues<dim> fe_values(fe,
                          QGauss<dim>(1),
                          update_inverse_jacobians | update_quadrature_points);

  for (const auto &cell : dof.active_cell_iterators())
    {
      fe_values.reinit(cell);
      deallog << "Data on cell " << cell->id()
              << ": p=" << fe_values.quadrature_point(0)
              << " jac=" << Tensor<2, dim>(fe_values.inverse_jacobian(0))
              << std::endl;
    }
  deallog << std::endl;

  return 0;
}
