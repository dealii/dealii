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


// Check precomputed mapping information against mapping
// information on the fly for a vector valued FE and MappingCartesian.
// This is a modified version of point_evaluation_13.cc.

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_cartesian.h>

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



template <int dim>
void
test()
{
  Triangulation<dim> tria;

  GridGenerator::subdivided_hyper_cube(tria, 2, 0, 1);

  MappingCartesian<dim> mapping;
  deallog << "Cartesian linear mapping" << std::endl;

  std::vector<Point<dim>> unit_points;
  for (unsigned int i = 0; i < 13; ++i)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = static_cast<double>(i) / 17. + 0.015625 * d;
      unit_points.push_back(p);
    }

  FESystem<dim>   fe(FE_Q<dim>(2), dim);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  Vector<double> vector(dof_handler.n_dofs());

  // mapping info object for precomputed mapping
  NonMatching::MappingInfo<dim, dim> mapping_info(mapping,
                                                  update_values |
                                                    update_gradients);
  // The first evaluator will use the precomputed mapping
  FEPointEvaluation<dim, dim> evaluator1(mapping_info, fe);

  // The second evaluator will compute the mapping on the fly
  FEPointEvaluation<dim, dim> evaluator2(mapping,
                                         fe,
                                         update_values | update_gradients);

  VectorTools::interpolate(mapping, dof_handler, MyFunction<dim>(), vector);

  std::vector<double> solution_values(fe.dofs_per_cell);
  std::vector<double> solution_values2(fe.dofs_per_cell);

  const FEValuesExtractors::Vector extractor(0);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell->get_dof_values(vector,
                           solution_values.begin(),
                           solution_values.end());

      mapping_info.reinit(cell, unit_points);
      evaluator1.evaluate(solution_values,
                          EvaluationFlags::values | EvaluationFlags::gradients);

      evaluator2.reinit(cell, unit_points);
      evaluator2.evaluate(solution_values,
                          EvaluationFlags::values | EvaluationFlags::gradients);

      deallog << "Cell with center " << cell->center(true) << std::endl;
      for (unsigned int i = 0; i < unit_points.size(); ++i)
        deallog
          << mapping.transform_unit_to_real_cell(cell, unit_points[i]) << ": "
          << evaluator1.get_value(i) << " error value "
          << (evaluator1.get_value(i) - evaluator2.get_value(i)).norm()
          << " error grad "
          << (evaluator1.get_gradient(i) - evaluator2.get_gradient(i)).norm()
          << std::endl;
      deallog << std::endl;

      for (unsigned int i = 0; i < unit_points.size(); ++i)
        {
          evaluator1.submit_value(evaluator1.get_value(i), i);
          evaluator1.submit_gradient(evaluator1.get_gradient(i), i);

          evaluator2.submit_value(evaluator2.get_value(i), i);
          evaluator2.submit_gradient(evaluator2.get_gradient(i), i);
        }

      evaluator1.test_and_sum(solution_values,
                              EvaluationFlags::values |
                                EvaluationFlags::gradients);
      evaluator2.test_and_sum(solution_values2,
                              EvaluationFlags::values |
                                EvaluationFlags::gradients);

      for (unsigned int i = 0; i < solution_values.size(); ++i)
        deallog << solution_values[i] - solution_values2[i] << ' ';
      deallog << std::endl;
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  test<2>();
  test<3>();
}
