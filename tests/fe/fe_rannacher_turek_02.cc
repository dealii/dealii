// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Interfaces being tested
#include <deal.II/base/polynomials_rannacher_turek.h>

#include <deal.II/fe/fe_rannacher_turek.h>
// Interfaces needed for testing
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


// Regression test for some values of Rannacher Turek element.

void
test_values()
{
  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria, 0.0, 1.0);

  FE_RannacherTurek<2> fe;
  DoFHandler<2>        dofh(tria);
  dofh.distribute_dofs(fe);

  QGauss<2>   quadrature(8);
  FEValues<2> fev(fe,
                  quadrature,
                  update_values | update_gradients | update_JxW_values);
  fev.reinit(dofh.begin_active());

  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
    {
      for (unsigned int q = 0; q < quadrature.size(); ++q)
        {
          deallog << fev.shape_value(i, q) << ' ';
          for (unsigned int d = 0; d < 2; ++d)
            {
              deallog << fev.shape_grad(i, q)[d] << ' ';
            }
        }
      deallog << std::endl;
    }
}

int
main()
{
  initlog();

  test_values();

  return 0;
}
