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



// Test MatrixFreeTools::compute_matrix():
//   (1) compute individual components
//   (2) compute with multiple FEEvaluation instances
//

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include <vector>

#include "../tests.h"

using namespace dealii;

template <unsigned int n_components,
          int          dim,
          typename Number,
          typename VectorizedArrayType>
void
run(const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const unsigned int first_selected_component)
{
  const unsigned int n_dofs = matrix_free.get_dof_handler().n_dofs();
  FullMatrix<double> matrix(n_dofs, n_dofs);

  FullMatrix<double> scaling(3, 3);
  for (unsigned int i = 0, c = 1; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j, ++c)
      scaling[i][j] = c;

  MatrixFreeTools::
    compute_matrix<dim, -1, 0, n_components, double, VectorizedArray<double>>(
      matrix_free,
      matrix_free.get_affine_constraints(),
      matrix,
      [&](auto &phi) {
        phi.evaluate(EvaluationFlags::values);
        for (const auto q : phi.quadrature_point_indices())
          {
            auto value = phi.get_value(q);

            auto result = value;
            result      = 0.0;

            for (unsigned int i = 0; i < n_components; ++i)
              for (unsigned int j = 0; j < n_components; ++j)
                if constexpr (n_components > 1)
                  result[i] += scaling[i + first_selected_component]
                                      [j + first_selected_component] *
                               value[j];
                else
                  result += scaling[i + first_selected_component]
                                   [j + first_selected_component] *
                            value;

            phi.submit_value(result, q);
          }
        phi.integrate(EvaluationFlags::values);
      },
      0,
      0,
      first_selected_component);

  matrix.print_formatted(deallog.get_file_stream(), 5, false, 10, "0.00000");
  deallog << std::endl;
}



template <int dim>
void
test_0()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FESystem<dim>(FE_Q<dim>(2), 2, FE_Q<dim>(1), 1));
  DoFRenumbering::component_wise(dof_handler);

  QGauss<dim> quadrature(2);

  MappingQ1<dim> mapping;

  AffineConstraints<double> constraints;

  MatrixFree<dim, double> matrix_free;

  matrix_free.reinit(mapping, dof_handler, constraints, quadrature);

  run<1>(matrix_free, 0u);
  run<1>(matrix_free, 1u);
  run<1>(matrix_free, 2u);

  run<2>(matrix_free, 0u);
}



template <int dim>
void
test_1()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FESystem<dim>(FE_Q<dim>(2), 2));
  DoFRenumbering::component_wise(dof_handler);

  QGauss<dim> quadrature(2);

  MappingQ1<dim> mapping;

  AffineConstraints<double> constraints;

  MatrixFree<dim, double> matrix_free;

  matrix_free.reinit(mapping, dof_handler, constraints, quadrature);

  MatrixFreeTools::internal::
    ComputeMatrixScratchData<dim, VectorizedArray<double>, false>
      data_cell;

  data_cell.dof_numbers               = {0, 0};
  data_cell.quad_numbers              = {0, 0};
  data_cell.n_components              = {1, 1};
  data_cell.first_selected_components = {0, 1};
  data_cell.batch_type                = {0, 0};

  data_cell.op_create =
    [&](const std::pair<unsigned int, unsigned int> &range) {
      std::vector<
        std::unique_ptr<FEEvaluationData<dim, VectorizedArray<double>, false>>>
        phi;

      phi.emplace_back(std::make_unique<FEEvaluation<dim, -1, 0, 1, double>>(
        matrix_free, range, 0, 0, 0));

      phi.emplace_back(std::make_unique<FEEvaluation<dim, -1, 0, 1, double>>(
        matrix_free, range, 0, 0, 1));

      return phi;
    };

  data_cell.op_reinit = [](auto &phi, const unsigned batch) {
    static_cast<FEEvaluation<dim, -1, 0, 1, double> &>(*phi[0]).reinit(batch);
    static_cast<FEEvaluation<dim, -1, 0, 1, double> &>(*phi[1]).reinit(batch);
  };

  data_cell.op_compute = [&](auto &phi) {
    auto &phi_0 = static_cast<FEEvaluation<dim, -1, 0, 1, double> &>(*phi[0]);
    auto &phi_1 = static_cast<FEEvaluation<dim, -1, 0, 1, double> &>(*phi[1]);

    phi_0.evaluate(EvaluationFlags::values);
    phi_1.evaluate(EvaluationFlags::values);
    for (const auto q : phi_0.quadrature_point_indices())
      {
        auto value_0 = phi_0.get_value(q);
        auto value_1 = phi_1.get_value(q);

        phi_0.submit_value(value_0 * 1.0 + value_1 * 2.0, q);
        phi_1.submit_value(value_0 * 4.0 + value_1 * 5.0, q);
      }
    phi_0.integrate(EvaluationFlags::values);
    phi_1.integrate(EvaluationFlags::values);
  };

  const unsigned int n_dofs = dof_handler.n_dofs();
  FullMatrix<double> matrix(n_dofs, n_dofs);

  MatrixFreeTools::internal::compute_matrix(
    matrix_free, constraints, data_cell, {}, {}, matrix);

  matrix.print_formatted(deallog.get_file_stream(), 5, false, 10, "0.00000");
  deallog << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  test_0<1>(); // individual components
  test_1<1>(); // multiple FEEvaluation instances
}
