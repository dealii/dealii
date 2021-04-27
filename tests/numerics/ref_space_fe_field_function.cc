// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2021 by the deal.II authors
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


/*
 * Test the class RefSpaceFEFieldFunction,
 * by setting up a single cell triangulation, a DoFHandler, and a vector with
 * solution values. Create an RefSpaceFEFieldFunction object, call the functions
 * value, gradient, and hessian, and print their values to deallog.
 */

#include <deal.II/base/point.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/fe_field_function.h>

#include "../tests.h"


using namespace dealii;

template <int dim>
void
run_test()
{
  deallog << "dim = " << dim << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);

  const FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  const unsigned int n_dofs = dof_handler.n_dofs();

  // We want to check the value/gradient/hessian at the center of the cell.
  const typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active();
  const Point<dim> cell_center = cell->center();

  // Choose some solution which have mostly non-zero values for
  // value/gradient/hessian.
  Vector<double> solution;
  solution.reinit(n_dofs);
  for (unsigned int i = 0; i < n_dofs; i++)
    solution[i] = std::pow(i + 1, 2);

  Functions::RefSpaceFEFieldFunction<dim> function(dof_handler, solution);
  function.set_active_cell(cell);

  deallog << "value = " << function.value(cell_center) << std::endl;
  deallog << "gradient = " << function.gradient(cell_center) << std::endl;
  deallog << "hessian = " << function.hessian(cell_center) << std::endl;

  deallog << std::endl;
}



int
main()
{
  initlog();
  run_test<1>();
  run_test<2>();
  run_test<3>();
}
