// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check FEPointEvaluation for dim != spacedim, otherwise same test as
// point_evaluation_03

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"



template <int dim>
class MyFunction : public Function<dim>
{
public:
  MyFunction()
    : Function<dim>(dim)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const override
  {
    return component + p[component];
  }
};



template <int dim, int spacedim>
void
test(const unsigned int degree)
{
  Triangulation<dim, spacedim> tria;

  if constexpr (dim + 1 == spacedim)
    GridGenerator::hyper_sphere(tria, Point<spacedim>(), 1.);
  else if constexpr (dim == 1)
    GridGenerator::hyper_rectangle(tria, Point<dim>(), Point<dim>(1.2));
  else
    static_assert("Case not implemented");

  MappingQ<dim, spacedim> mapping(degree);
  deallog << "Mapping of degree " << degree << std::endl;

  std::vector<Point<dim>> unit_points;
  for (unsigned int i = 0; i < 13; ++i)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = static_cast<double>(i) / 17. + 0.015625 * d;
      unit_points.push_back(p);
    }

  FESystem<dim, spacedim> fe(FE_Q<dim, spacedim>(degree), spacedim);
  FEValues<dim, spacedim> fe_values(mapping,
                                    fe,
                                    Quadrature<dim>(unit_points),
                                    update_values | update_gradients);

  DoFHandler<dim, spacedim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  Vector<double> vector(dof_handler.n_dofs());

  FEPointEvaluation<spacedim, dim, spacedim> evaluator(
    mapping, fe, update_values | update_gradients);

  VectorTools::interpolate(mapping,
                           dof_handler,
                           MyFunction<spacedim>(),
                           vector);

  std::vector<double>              solution_values(fe.dofs_per_cell);
  std::vector<Tensor<1, spacedim>> function_values(unit_points.size());
  std::vector<Tensor<2, spacedim>> function_gradients(unit_points.size());

  const FEValuesExtractors::Vector extractor(0);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      fe_values[extractor].get_function_values(vector, function_values);
      fe_values[extractor].get_function_gradients(vector, function_gradients);

      cell->get_dof_values(vector,
                           solution_values.begin(),
                           solution_values.end());

      evaluator.reinit(cell, unit_points);
      evaluator.evaluate(solution_values,
                         EvaluationFlags::values | EvaluationFlags::gradients);

      deallog << "Cell with center " << cell->center(true) << std::endl;
      for (unsigned int i = 0; i < function_values.size(); ++i)
        deallog << mapping.transform_unit_to_real_cell(cell, unit_points[i])
                << ": " << evaluator.get_value(i) << " error value "
                << (function_values[i] - evaluator.get_value(i)).norm()
                << " error grad "
                << (evaluator.get_gradient(i) - function_gradients[i]).norm()
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
        deallog << i << ' ';
      deallog << std::endl;
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  test<1, 2>(2);
  test<1, 3>(6);
  test<2, 3>(5);
}
