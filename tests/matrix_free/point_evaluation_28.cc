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


// check FEPointEvaluation integrate_add()

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
  for (unsigned int i = 0; i < Utilities::pow(2, dim); ++i)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = static_cast<Number>(i) / 17. + 0.015625 * d;
      unit_points.push_back(p);
    }

  QGauss<dim> q_gauss(2);

  Quadrature<dim> quad(unit_points, q_gauss.get_weights());

  FE_Q<dim> fe(degree);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  Vector<Number> vector(dof_handler.n_dofs());

  NonMatching::MappingInfo mapping_info(mapping,
                                        update_values | update_gradients |
                                          update_JxW_values);

  std::vector<Quadrature<dim>> quad_vec_cell;
  quad_vec_cell.reserve(tria.n_cells());

  for (const auto &cell : tria.active_cell_iterators())
    quad_vec_cell.push_back(quad);

  mapping_info.reinit_cells(tria.active_cell_iterators(), quad_vec_cell);

  FEPointEvaluation<1, dim, dim, Number> evaluator(mapping_info, fe);

  Tensor<1, dim, Number> exponents;
  exponents[0] = 1.;
  VectorTools::interpolate(mapping,
                           dof_handler,
                           Functions::Monomial<dim, Number>(exponents),
                           vector);

  std::vector<Number> solution_values(fe.dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell->get_dof_values(vector,
                           solution_values.begin(),
                           solution_values.end());

      evaluator.reinit(cell->active_cell_index());
      evaluator.evaluate(solution_values,
                         EvaluationFlags::values | EvaluationFlags::gradients);

      for (const unsigned int i : evaluator.quadrature_point_indices())
        {
          evaluator.submit_value(evaluator.get_value(i), i);
          evaluator.submit_gradient(evaluator.get_gradient(i), i);
        }

      evaluator.integrate(solution_values,
                          EvaluationFlags::values | EvaluationFlags::gradients);

      deallog << "Cell with center " << cell->center(true) << std::endl;
      for (const auto i : solution_values)
        deallog << i << ' ';
      deallog << std::endl;

      evaluator.integrate(solution_values,
                          EvaluationFlags::values | EvaluationFlags::gradients,
                          true);

      for (const auto i : solution_values)
        deallog << i << ' ';
      deallog << std::endl;
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  test<2>(1);
  test<2>(2);
}
