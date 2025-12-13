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


// check FEPointEvaluation for a case where the evaluator sets two component,
// but the initialization is with a scalar FE_DGQ and MappingQ

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <iostream>

#include "../tests.h"



template <int dim, typename Number = double>
void
test(const FiniteElement<dim> &fe)
{
  deallog << "Element " << fe.get_name() << std::endl;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1.2, 1.2);

  MappingQ<dim> mapping(fe.degree);

  std::vector<Point<dim>> unit_points;
  for (unsigned int i = 0; i < 13; ++i)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = static_cast<Number>(i) / 17. + 0.015625 * d;
      unit_points.push_back(p);
    }

  std::vector<Number> vector(fe.dofs_per_cell * dim);
  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
      vector[i + d * fe.dofs_per_cell] = (0.5 + 0.2 * d) * i + d;

  FEPointEvaluation<dim, dim, dim, Number> evaluator(
    mapping, fe, update_values | update_gradients);

  const auto cell = tria.begin_active();
  evaluator.reinit(cell, unit_points);
  evaluator.evaluate(vector,
                     EvaluationFlags::values | EvaluationFlags::gradients);

  for (unsigned int i = 0; i < unit_points.size(); ++i)
    deallog << "At "
            << mapping.transform_unit_to_real_cell(cell, unit_points[i])
            << "  value " << evaluator.get_value(i) << "  gradient "
            << evaluator.get_gradient(i) << std::endl;
  deallog << std::endl;

  for (unsigned int i = 0; i < unit_points.size(); ++i)
    {
      evaluator.submit_value(evaluator.get_value(i), i);
      evaluator.submit_gradient(evaluator.get_gradient(i), i);
    }

  evaluator.test_and_sum(vector,
                         EvaluationFlags::values | EvaluationFlags::gradients);

  deallog << "Testing:" << std::endl;
  for (const Number i : vector)
    deallog << i << ' ';
  deallog << std::endl << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  test(FE_DGQ<2>(2));
  test(FE_DGQ<3>(2));
  test(FE_Q<3>(2));
}
