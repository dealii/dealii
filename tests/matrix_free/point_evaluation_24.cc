// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2023 by the deal.II authors
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


// Check face path of FEPointEvaluation with precomputed mapping
// information for a vector valued FE and MappingQ.


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

  std::vector<Point<dim - 1>> unit_points;
  for (unsigned int i = 0; i < (dim > 1 ? 13 : 1); ++i)
    {
      Point<dim - 1> p;
      for (unsigned int d = 0; d < dim - 1; ++d)
        p[d] = static_cast<double>(i) / 17. + 0.015625 * d;
      unit_points.push_back(p);
    }

  FESystem<dim>   fe(FE_Q<dim>(2), dim);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  Vector<double> vector(dof_handler.n_dofs());

  NonMatching::MappingInfo<dim, dim> mapping_info(mapping,
                                                  update_values |
                                                    update_gradients);

  std::vector<std::vector<Quadrature<dim - 1>>> quad_vec_faces(
    tria.n_active_cells());

  for (const auto &cell : tria.active_cell_iterators())
    {
      for (const auto f : cell->face_indices())
        {
          quad_vec_faces[cell->active_cell_index()].push_back(
            Quadrature<dim - 1>(unit_points));
        }
    }

  mapping_info.reinit_faces(tria.active_cell_iterators(), quad_vec_faces);

  FEPointEvaluation<dim, dim> evaluator(mapping_info, fe);

  VectorTools::interpolate(mapping, dof_handler, MyFunction<dim>(), vector);

  std::vector<double> solution_values_in(fe.dofs_per_cell);
  std::vector<double> solution_values_out(fe.dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell->get_dof_values(vector,
                           solution_values_in.begin(),
                           solution_values_in.end());

      deallog << "Cell with center " << cell->center(true) << std::endl;
      for (const auto f : cell->face_indices())
        {
          evaluator.reinit(cell->active_cell_index(), f);

          evaluator.evaluate(solution_values_in,
                             EvaluationFlags::values |
                               EvaluationFlags::gradients);
          deallog << "face number " << f << std::endl;
          for (unsigned int i = 0; i < unit_points.size(); ++i)
            deallog << "value " << evaluator.get_value(i) << " gradient "
                    << evaluator.get_gradient(i) << std::endl;
          deallog << std::endl;

          for (unsigned int i = 0; i < unit_points.size(); ++i)
            {
              evaluator.submit_value(evaluator.get_value(i), i);
              evaluator.submit_gradient(evaluator.get_gradient(i), i);
            }

          evaluator.test_and_sum(solution_values_out,
                                 EvaluationFlags::values |
                                   EvaluationFlags::gradients);

          for (double i : solution_values_out)
            deallog << i << ' ';
          deallog << std::endl;
        }
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  test<1>(1);
  test<1>(2);
  test<2>(1);
  test<2>(2);
  test<2>(6);
  test<3>(5);
}
