// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// turns out that cylinder_shell_01 wasn't enough: measure just takes the
// positive value it computes. we need to ask a FEValues object for the JxW
// values to get the sign of the Jacobian


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



int
main()
{
  initlog();
  deallog << std::setprecision(2);

  // generate a hyperball in 3d
  Triangulation<3> tria;
  GridGenerator::cylinder_shell(tria, 1, .8, 1);

  FE_Q<3>       fe(1);
  DoFHandler<3> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  QMidpoint<3> q;
  FEValues<3>  fe_values(fe, q, update_JxW_values);

  // make sure that all cells have positive
  // volume
  for (DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      fe_values.reinit(cell);
      deallog << cell << ' ' << fe_values.JxW(0) << std::endl;

      Assert(fe_values.JxW(0) > 0, ExcInternalError());
    }
}
