// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check FEPointEvaluation::quadrature_point(), FEPointEvaluation::unit_point(),
// FEPointEvaluation::jacobian(), FEPointEvaluation::inverse_jacobian(),
// FEPointEvaluation::get_unit_gradient().

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"



template <int dim>
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
  for (unsigned int i = 0; i < 7; ++i)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = static_cast<double>(i) / 17. + 0.015625 * d;
      unit_points.push_back(p);
    }

  FE_Q<dim> fe(degree);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  Vector<double> vector(dof_handler.n_dofs());

  Tensor<1, dim> exponents;
  exponents[0] = 1.;
  VectorTools::interpolate(mapping,
                           dof_handler,
                           Functions::Monomial<dim>(exponents),
                           vector);

  FEPointEvaluation<1, dim> evaluator(mapping,
                                      fe,
                                      update_values | update_gradients |
                                        update_quadrature_points |
                                        update_jacobians);
  std::vector<double>       solution_values(fe.dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell->get_dof_values(vector,
                           solution_values.begin(),
                           solution_values.end());

      evaluator.reinit(cell, unit_points);
      evaluator.evaluate(solution_values,
                         EvaluationFlags::values | EvaluationFlags::gradients);

      deallog << "Cell with center " << cell->center(true) << std::endl;
      for (unsigned int i = 0; i < unit_points.size(); ++i)
        deallog << "unit point " << unit_points[i] << std::endl
                << "unit point via evaluator: " << evaluator.unit_point(i)
                << std::endl
                << "real point: " << evaluator.quadrature_point(i) << std::endl
                << "jacobian: " << Tensor<2, dim>(evaluator.jacobian(i))
                << std::endl
                << "inverse jacobian: "
                << Tensor<2, dim>(evaluator.inverse_jacobian(i)) << std::endl
                << std::endl;
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  test<1>(3);
  test<2>(2);
  test<2>(6);
  test<3>(5);
}
