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


// Check face path of FEPointEvaluation with dim < spacedim, otherwise similar
// to point_evaluation_23


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
    : Function<dim>()
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

  std::vector<Point<dim - 1>> unit_points;
  for (unsigned int i = 0; i < (dim > 1 ? 13 : 1); ++i)
    {
      Point<dim - 1> p;
      for (unsigned int d = 0; d < dim - 1; ++d)
        p[d] = static_cast<double>(i) / 17. + 0.015625 * d;
      unit_points.push_back(p);
    }

  FE_Q<dim, spacedim>       fe(degree);
  DoFHandler<dim, spacedim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  Vector<double> vector(dof_handler.n_dofs());

  NonMatching::MappingInfo<dim, spacedim> mapping_info(mapping,
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

  FEFacePointEvaluation<1, dim, spacedim> evaluator(mapping_info, fe);

  VectorTools::interpolate(mapping,
                           dof_handler,
                           MyFunction<spacedim>(),
                           vector);

  FEFaceValues<dim, spacedim> fe_val(mapping,
                                     fe,
                                     Quadrature<dim - 1>(unit_points),
                                     update_values | update_gradients);

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
          fe_val.reinit(cell, f);

          evaluator.evaluate(solution_values_in,
                             EvaluationFlags::values |
                               EvaluationFlags::gradients);
          deallog << "face number " << f << std::endl;
          for (unsigned int i = 0; i < unit_points.size(); ++i)
            deallog << "value " << evaluator.get_value(i) << " gradient "
                    << evaluator.get_gradient(i) << std::endl;
          deallog << std::endl;
          deallog << "reference" << std::endl;
          for (unsigned int i = 0; i < unit_points.size(); ++i)
            {
              double              val = 0;
              Tensor<1, spacedim> grad;
              for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
                {
                  val += solution_values_in[j] * fe_val.shape_value(j, i);
                  grad += solution_values_in[j] * fe_val.shape_grad(j, i);
                }
              deallog << "value " << val << " gradient " << grad << std::endl;
            }
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

  test<1, 2>(2);
  test<1, 3>(2);
  test<2, 3>(3);
}
