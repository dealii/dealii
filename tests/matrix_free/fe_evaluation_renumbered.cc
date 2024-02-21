// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


// Test FEEvaluation::reinit() that thaes std::array of indices.



template <int dim>
class RightHandSideFunction : public Function<dim>
{
public:
  RightHandSideFunction(const unsigned int n_components = 1)
    : Function<dim>(n_components)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    return p[component % dim] * p[component % dim];
  }
};

template <int dim, int fe_degree, int n_points, typename Number>
void
test(const unsigned int n_refinements, const unsigned int geometry_type)
{
  using VectorizedArrayType = VectorizedArray<Number>;

  using VectorType = Vector<Number>;

  Triangulation<dim> tria;

  if (geometry_type == 0)
    GridGenerator::hyper_cube(tria);
  else if (geometry_type == 1)
    GridGenerator::hyper_ball(tria);
  else
    AssertThrow(false, ExcNotImplemented());

  tria.refine_global(n_refinements);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  MappingQ<dim> mapping(1);
  QGauss<1>     quad(n_points);

  AffineConstraints<Number> constraint;

  using MF = MatrixFree<dim, Number, VectorizedArrayType>;

  typename MF::AdditionalData additional_data;
  additional_data.mapping_update_flags =
    update_values | update_quadrature_points;
  additional_data.tasks_parallel_scheme = MF::AdditionalData::none;

  MF matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraint, quad, additional_data);

  VectorType src, dst;

  matrix_free.initialize_dof_vector(src);
  matrix_free.initialize_dof_vector(dst);

  VectorTools::interpolate(mapping,
                           dof_handler,
                           RightHandSideFunction<dim>(1),
                           src);

  std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> indices;

  for (unsigned int cell = 0; cell < matrix_free.n_cell_batches(); ++cell)
    {
      for (unsigned int i = 0;
           i < matrix_free.n_active_entries_per_cell_batch(cell);
           ++i)
        indices.emplace_back(indices.size(), cell, i);
    }

  // if(false)
  std::sort(indices.begin(), indices.end(), [](const auto &a, const auto &b) {
    if (std::get<2>(a) != std::get<2>(b))
      return std::get<2>(a) < std::get<2>(b);

    return std::get<1>(a) < std::get<1>(b);
  });

  // collect reference data with normal reinit function

  // quadrature points -> mapping
  std::vector<std::vector<Point<dim>>> quadrature_points_ref;

  // dof values -> read/write operation
  std::vector<std::vector<Number>> values_ref;

  matrix_free.template cell_loop<VectorType, VectorType>(
    [&](const auto &matrix_free, auto &dst, const auto &src, auto range) {
      FEEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>
        phi(matrix_free);

      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          phi.reinit(cell);
          phi.read_dof_values(src);

          for (unsigned int v = 0;
               v < matrix_free.n_active_entries_per_cell_batch(cell);
               ++v)
            {
              std::vector<Point<dim>> points;

              for (unsigned int i = 0; i < phi.n_q_points; ++i)
                {
                  auto temp_v = phi.quadrature_point(i);

                  Point<dim> temp;
                  for (unsigned int d = 0; d < dim; ++d)
                    temp[d] = temp_v[d][v];

                  points.emplace_back(temp);
                }

              quadrature_points_ref.emplace_back(points);

              std::vector<Number> values;

              for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
                values.emplace_back(phi.get_dof_value(i)[v]);
              values_ref.emplace_back(values);
            }
        }
    },
    dst,
    src);

  // test variable reinit() function
  {
    FEEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType> phi(
      matrix_free);

    for (unsigned int v = 0; v < indices.size();
         v += VectorizedArrayType::size())
      {
        std::array<unsigned int, VectorizedArrayType::size()> indices_;

        indices_.fill(numbers::invalid_unsigned_int);

        const unsigned int n_lanes_filled =
          std::min(v + VectorizedArrayType::size(), indices.size()) - v;

        for (unsigned int i = v, c = 0; i < v + n_lanes_filled; ++i, ++c)
          indices_[c] = std::get<1>(indices[i]) * VectorizedArrayType::size() +
                        std::get<2>(indices[i]);

        phi.reinit(indices_);
        phi.read_dof_values(src);

        for (unsigned int i = v, c = 0; i < v + n_lanes_filled; ++i, ++c)
          {
            std::vector<Point<dim>> points;

            for (unsigned int i = 0; i < phi.n_q_points; ++i)
              {
                auto temp_v = phi.quadrature_point(i);

                Point<dim> temp;
                for (unsigned int d = 0; d < dim; ++d)
                  temp[d] = temp_v[d][c];

                points.emplace_back(temp);
              }

            // perform comparison

            Assert(points == quadrature_points_ref[std::get<0>(indices[i])],
                   ExcInternalError());

            std::vector<Number> values;

            for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
              values.emplace_back(phi.get_dof_value(i)[c]);

            // perform comparison
            Assert(values == values_ref[std::get<0>(indices[i])],
                   ExcInternalError());
          }
      }
  }

  deallog << "OK!" << std::endl;
}

int
main()
{
  initlog();
  test<2, 1, 2, double>(3, 0);
  test<2, 1, 2, double>(3, 1);
}
