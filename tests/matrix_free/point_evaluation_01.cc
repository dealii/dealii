// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check FEPointEvaluation for scalar FE_Q and MappingQ by comparing to
// the output of FEValues with the same settings

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"



template <int dim, typename Number = double>
void
test(const unsigned int degree)
{
  Triangulation<dim> tria;

  if (dim > 1)
    GridGenerator::hyper_shell(tria, Point<dim>(), 0.5, 1, 6);
  else
    GridGenerator::subdivided_hyper_cube(tria, 2, 0, 1);

  MappingQ<dim> mapping(degree);
  deallog << "Mapping of degree " << degree << std::endl;

  std::vector<Point<dim>> unit_points;
  for (unsigned int i = 0; i < 13; ++i)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = static_cast<Number>(i) / 17. + 0.015625 * d;
      unit_points.push_back(p);
    }

  FE_Q<dim>     fe(degree);
  FEValues<dim> fe_values(mapping,
                          fe,
                          Quadrature<dim>(unit_points),
                          update_values | update_gradients);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  Vector<Number> vector(dof_handler.n_dofs());

  FEPointEvaluation<1, dim, dim, Number> evaluator(
    mapping, fe, update_values | update_gradients);

  Tensor<1, dim, Number> exponents;
  exponents[0] = 1.;
  VectorTools::interpolate(mapping,
                           dof_handler,
                           Functions::Monomial<dim, Number>(exponents),
                           vector);

  std::vector<Number>                 solution_values(fe.dofs_per_cell);
  std::vector<Number>                 function_values(unit_points.size());
  std::vector<Tensor<1, dim, Number>> function_gradients(unit_points.size());

  // For float numbers that are sensitive to roundoff in the numdiff
  // tolerances (absolute 1e-8), we multiply by 1e-3 to ensure that the test
  // remains robust
  const double factor_float = std::is_same_v<Number, float> ? 0.001 : 1.;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      fe_values.get_function_values(vector, function_values);
      fe_values.get_function_gradients(vector, function_gradients);

      cell->get_dof_values(vector,
                           solution_values.begin(),
                           solution_values.end());

      evaluator.reinit(cell, unit_points);
      evaluator.evaluate(solution_values,
                         EvaluationFlags::values | EvaluationFlags::gradients);

      deallog << "Cell with center " << cell->center(true) << std::endl;
      for (unsigned int i = 0; i < function_values.size(); ++i)
        deallog << mapping.transform_unit_to_real_cell(cell, unit_points[i])
                << ": " << factor_float * evaluator.get_value(i)
                << " error value "
                << factor_float * (function_values[i] - evaluator.get_value(i))
                << " error grad "
                << factor_float *
                     (evaluator.get_gradient(i) - function_gradients[i]).norm()
                << std::endl;
      deallog << std::endl;

      for (unsigned int i = 0; i < unit_points.size(); ++i)
        {
          evaluator.submit_value(evaluator.get_value(i), i);
          evaluator.submit_gradient(evaluator.get_gradient(i), i);
        }

      evaluator.test_and_sum(solution_values,
                             EvaluationFlags::values |
                               EvaluationFlags::gradients);

      for (const auto i : solution_values)
        deallog << factor_float * i << ' ';
      deallog << std::endl;
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  test<1>(3);
  test<2>(1);
  test<2>(2);
  test<2>(6);
  test<3>(1);
  test<3>(5);

  test<3, float>(5);
}
