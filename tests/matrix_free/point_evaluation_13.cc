// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


// Check FEPointEvaluation for vector-valued FE_Q and MappingCartesian by
// comparing to the output of FEValues with the same settings.
// This is a modified version of point_evaluation_03.cc.

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
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
  using namespace dealii;
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

  FESystem<dim> fe(FE_Q<dim>(2), dim);
  FEValues<dim> fe_values(mapping,
                          fe,
                          Quadrature<dim>(unit_points),
                          update_values | update_gradients);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  Vector<double> vector(dof_handler.n_dofs());

  FEPointEvaluation<dim, dim> evaluator(mapping,
                                        fe,
                                        update_values | update_gradients);

  VectorTools::interpolate(mapping, dof_handler, MyFunction<dim>(), vector);

  std::vector<double>         solution_values(fe.dofs_per_cell);
  std::vector<Tensor<1, dim>> function_values(unit_points.size());
  std::vector<Tensor<2, dim>> function_gradients(unit_points.size());

  FEValuesExtractors::Vector extractor(0);

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

      evaluator.integrate(solution_values,
                          EvaluationFlags::values | EvaluationFlags::gradients);

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

  test<2>();
  test<3>();
}
