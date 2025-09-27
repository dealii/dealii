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



// Test ShapeData for FE_SimplexP and QGaussSimplex

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
class ExactSolution : public Function<dim>
{
public:
  ExactSolution()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int /*component*/ = 0) const
  {
    return p[0] + 2 * p[1];
  }
};

template <int dim, typename Number>
class FEEvaluationDummy
{
public:
  FEEvaluationDummy(
    const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info)
    : shape_info(shape_info)
  {
    values.resize(shape_info.n_q_points);
    gradients.resize(shape_info.n_q_points * dim);
  }

  void
  evaluate(const Vector<Number> &dof_values)
  {
    const unsigned int n_dofs     = shape_info.dofs_per_component_on_cell;
    const unsigned int n_q_points = shape_info.n_q_points;

    for (unsigned int q = 0, c = 0; q < shape_info.n_q_points; ++q)
      {
        values[q] = 0.0;
        for (unsigned int i = 0; i < n_dofs; ++i, ++c)
          values[q] += shape_info.data[0].shape_values[c] * dof_values[i];
      }

    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          gradients[n_q_points * d + q] = 0.0;
          for (unsigned int i = 0; i < n_dofs; ++i)
            gradients[n_q_points * d + q] +=
              shape_info.data[0]
                .shape_gradients[i * dim * n_q_points + q * dim + d] *
              dof_values[i];
        }
  }

  Number
  get_value(const unsigned int q_point) const
  {
    return values[q_point];
  }

  Tensor<1, dim, Number>
  get_gradient(const unsigned int q_point) const
  {
    Tensor<1, dim, Number> out;

    for (unsigned int c = 0; c < dim; ++c)
      out[c] = gradients[shape_info.n_q_points * c + q_point];

    return out;
  }

private:
  const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info;

  AlignedVector<Number> values;
  AlignedVector<Number> gradients;
};

template <int dim, int spacedim = dim, typename Number = double>
void
test(const FiniteElement<dim, spacedim> &fe)
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  MappingFE<dim> mapping(FE_SimplexP<dim>(1));

  QGaussSimplex<dim> quadrature(1);

  internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info(quadrature, fe);

  FEValues<dim> fe_values(mapping,
                          fe,
                          quadrature,
                          update_values | update_gradients |
                            update_inverse_jacobians);

  Vector<Number> src(dof_handler.n_dofs());
  VectorTools::interpolate(mapping, dof_handler, ExactSolution<dim>(), src);

  FEEvaluationDummy<dim, Number> fe_eval(shape_info);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);

      Vector<Number> dof_values(fe_values.dofs_per_cell);
      cell->get_dof_values(src, dof_values);
      fe_eval.evaluate(dof_values);

      // check values
      {
        std::vector<Number> values(fe_values.n_quadrature_points);
        fe_values.get_function_values(src, values);

        for (const auto q : fe_values.quadrature_point_indices())
          {
            Assert(std::abs(values[q] - fe_eval.get_value(q)) < 1e-10,
                   ExcMessage("Entries do not match!"));
          }
      }

      // check gradients
      {
        std::vector<Tensor<1, dim, Number>> gradients(
          fe_values.n_quadrature_points);
        fe_values.get_function_gradients(src, gradients);

        for (const auto q : fe_values.quadrature_point_indices())
          {
            // apply the transposed inverse Jacobian of the mapping
            Tensor<1, dim, Number> temp = fe_eval.get_gradient(q);
            Tensor<1, dim, Number> temp_transformed;
            for (int d = 0; d < dim; ++d)
              {
                Number sum = 0;
                for (int e = 0; e < dim; ++e)
                  sum += fe_values.inverse_jacobian(q)[e][d] * temp[e];
                temp_transformed[d] = sum;
              }

            for (int i = 0; i < dim; ++i)
              Assert(std::abs(gradients[q][i] - temp_transformed[i]) < 1e-10,
                     ExcMessage("Entries do not match!"));
          }
      }
    }
}

int
main()
{
  initlog();

  test<2>(FE_SimplexP<2>(2));
}
