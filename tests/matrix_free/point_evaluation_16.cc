// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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


// Check FEPointEvaluation::reinit with precomputed mapping
// information for a vector valued FE and MappingQ.
// This is a modified version of point_evaluation_03.cc.


#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
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



template <int dim>
void
test(const unsigned int degree)
{
  using namespace dealii;
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
        p[d] = static_cast<double>(i) / 17. + 0.015625 * d;
      unit_points.push_back(p);
    }

  FESystem<dim>   fe(FE_Q<dim>(degree), dim);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  Vector<double> vector(dof_handler.n_dofs());

  NonMatching::MappingInfo<dim, dim> mapping_info(mapping,
                                                  update_values |
                                                    update_gradients);

  FEPointEvaluation<dim, dim> evaluator(mapping_info, fe);
  FEPointEvaluation<dim, dim> evaluator2(mapping_info, fe);

  VectorTools::interpolate(mapping, dof_handler, MyFunction<dim>(), vector);

  std::vector<double> solution_values(fe.dofs_per_cell);
  std::vector<double> solution_values2(fe.dofs_per_cell);

  FEValuesExtractors::Vector extractor(0);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell->get_dof_values(vector,
                           solution_values.begin(),
                           solution_values.end());

      mapping_info.reinit(cell, unit_points);

      evaluator.evaluate(solution_values,
                         EvaluationFlags::values | EvaluationFlags::gradients);

      evaluator2.evaluate(solution_values,
                          EvaluationFlags::values | EvaluationFlags::gradients);

      deallog << "Cell with center " << cell->center(true) << std::endl;
      for (unsigned int i = 0; i < unit_points.size(); ++i)
        deallog
          << mapping.transform_unit_to_real_cell(cell, unit_points[i]) << ": "
          << evaluator.get_value(i) << " error value "
          << (evaluator.get_value(i) - evaluator2.get_value(i)).norm()
          << " error grad "
          << (evaluator.get_gradient(i) - evaluator2.get_gradient(i)).norm()
          << std::endl;
      deallog << std::endl;

      for (unsigned int i = 0; i < unit_points.size(); ++i)
        {
          evaluator.submit_value(evaluator.get_value(i), i);
          evaluator.submit_gradient(evaluator.get_gradient(i), i);

          evaluator2.submit_value(evaluator2.get_value(i), i);
          evaluator2.submit_gradient(evaluator2.get_gradient(i), i);
        }

      evaluator.integrate(solution_values,
                          EvaluationFlags::values | EvaluationFlags::gradients);
      evaluator2.integrate(solution_values2,
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

  test<2>(2);
  test<2>(6);
  test<3>(5);
}
