// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



/* Author: Markus Buerg, 2012 */
/* Purpose: Check FEFieldFunction for hp::DoFHandler. */



#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/fe_field_function.h>

#include "../tests.h"



template <int dim>
void
check()
{
  Triangulation<dim> triangulation;

  GridGenerator::subdivided_hyper_cube(triangulation, 2);
  hp::FECollection<dim> fe_collection;

  for (unsigned int i = 1; i <= triangulation.n_active_cells(); ++i)
    fe_collection.push_back(FE_Q<dim>(i));

  DoFHandler<dim> dof_handler(triangulation);

  dof_handler.distribute_dofs(fe_collection);

  Vector<double> vector(dof_handler.n_dofs());

  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    vector(i) = i;

  Functions::FEFieldFunction<dim> fe_field(dof_handler, vector);
  QGauss<dim>                     quadrature(5);

  deallog << "values:" << std::endl;

  std::vector<double> values(quadrature.size());

  fe_field.value_list(quadrature.get_points(), values);

  for (unsigned int q_point = 0; q_point < quadrature.size(); ++q_point)
    deallog << values[q_point] << std::endl;

  deallog << "gradients:" << std::endl;

  std::vector<Tensor<1, dim>> gradients(quadrature.size());

  fe_field.gradient_list(quadrature.get_points(), gradients);

  for (unsigned int q_point = 0; q_point < quadrature.size(); ++q_point)
    deallog << gradients[q_point] << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog << std::fixed;

  deallog.push("1d");
  check<1>();
  deallog.pop();
  deallog.push("2d");
  check<2>();
  deallog.pop();
  deallog.push("3d");
  check<3>();
  deallog.pop();
}
