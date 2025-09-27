// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test internal::store_affine_constraints() and get_affine_constraints()

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"


int
main()
{
  initlog();

  const unsigned int fe_degree = 1;
  constexpr int      dim       = 2;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0, 1, true);

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  QGauss<1> quadrature(fe_degree + 1);

  AffineConstraints<double> double_ac;
  double_ac.close();

  {
    typename MatrixFree<dim, double>::AdditionalData additional_data_double;
    MatrixFree<dim, double>                          matrix_free_double;
    matrix_free_double.reinit(MappingQ1<dim>{},
                              dof_handler,
                              double_ac,
                              quadrature,
                              additional_data_double);

    if (&matrix_free_double.get_affine_constraints() == &double_ac)
      deallog << "OK" << std::endl;
  }

  {
    typename MatrixFree<dim, float>::AdditionalData additional_data_float;
    MatrixFree<dim, float>                          matrix_free_float;
    matrix_free_float.reinit(MappingQ1<dim>{},
                             dof_handler,
                             double_ac,
                             quadrature,
                             additional_data_float);

    try
      {
        matrix_free_float.get_affine_constraints();
      }
    catch (...)
      {
        deallog << "OK" << std::endl;
      }
  }
}
