// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
//
//
// Fix bug with FEFaceValues initialization when used with FE_NedelecSZ

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec_sz.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <iostream>
#include <sstream>

#include "../tests.h"


template <int dim>
void
test(unsigned p)
{
  FE_NedelecSZ<dim>  fe(p);
  Triangulation<dim> triangulation;

  GridGenerator::hyper_cube(triangulation, -1., 1.);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  QGauss<dim - 1>   quadrature(3);
  FEFaceValues<dim> fe_face_values(fe,
                                   quadrature,
                                   update_values | update_gradients);
  fe_face_values.reinit(dof_handler.begin_active(), 0);

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  for (unsigned int p = 0; p < 3; ++p)
    {
      test<2>(p);
      test<3>(p);
    }
}
