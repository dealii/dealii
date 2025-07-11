// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check FEPointEvaluation lexicographic numbering

#include <deal.II/base/vectorization.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_update_flags.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <cstdlib>

#include "../tests.h"

using namespace dealii;

int
main(int argc, char *argv[])
{
  initlog();

  constexpr unsigned int dim        = 3;
  constexpr unsigned int degree     = 2;
  constexpr unsigned int components = 2;

  using Number                   = double;
  constexpr unsigned int n_lanes = VectorizedArray<Number>::size();

  const unsigned int cell_batch_index = 0;

  Triangulation<dim> triangulation;
  MappingQ<dim>      mapping(1);

  GridGenerator::subdivided_hyper_cube(triangulation, 4, -1, 1);

  auto dof_handler = DoFHandler<dim>(triangulation);

  FESystem<dim> fe_system(FE_Q<dim>(degree), components);

  dof_handler.distribute_dofs(fe_system);

  auto matrix_free = MatrixFree<dim, double, VectorizedArray<double>>();
  typename MatrixFree<dim, Number, VectorizedArray<Number>>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags =
    update_gradients | update_values | update_quadrature_points;
  additional_data.hold_all_faces_to_owned_cells = true;

  AffineConstraints<Number>   affine_constraints{};
  QGauss<dim>                 quadrature(degree + 1);
  std::vector<Tensor<1, dim>> normals(quadrature.size());
  for (auto &normal : normals)
    {
      const Number x = 1. / ((dim == 3) ? std::sqrt(3.) : std::sqrt(2.));
      normal[0]      = x;
      normal[1]      = x;
      if (dim == 3)
        normal[2] = x;
    }
  NonMatching::ImmersedSurfaceQuadrature<dim> surface_quadrature(
    quadrature.get_points(), quadrature.get_weights(), normals);

  // setup matrixfree object
  matrix_free.reinit(
    mapping, dof_handler, affine_constraints, quadrature, additional_data);

  NonMatching::MappingInfo<dim, dim, VectorizedArray<Number>> mapping_info(
    mapping,
    additional_data.mapping_update_flags | update_normal_vectors |
      update_JxW_values);

  std::vector<TriaActiveIterator<CellAccessor<dim>>> cell_vec(
    1, triangulation.begin_active());
  std::vector<NonMatching::ImmersedSurfaceQuadrature<dim>> quad_vec(
    1, surface_quadrature);
  mapping_info.reinit_surface(cell_vec, quad_vec);

  FEPointEvaluation<components, dim, dim, VectorizedArray<Number>>
    evaluator_surface(mapping_info, fe_system, 0, true);

  FEEvaluation<dim, -1, 0, components, Number, VectorizedArray<Number>>
    evaluator(matrix_free, 0, 0);

  LinearAlgebra::distributed::Vector<Number> src;
  matrix_free.initialize_dof_vector(src);
  for (unsigned int i = 0; i < src.size(); i++)
    src[i] = random_value<Number>();

  const auto dofs_per_cell = evaluator.dofs_per_cell;

  evaluator.reinit(cell_batch_index);
  evaluator.read_dof_values(src);

  for (unsigned int v = 0; v < 1; ++v)
    {
      evaluator_surface.reinit(0);

      StridedArrayView<const Number, n_lanes> strided_view(
        &evaluator.begin_dof_values()[0][v], dofs_per_cell);

      evaluator_surface.evaluate(StridedArrayView<const Number, n_lanes>(
                                   &evaluator.begin_dof_values()[0][v],
                                   dofs_per_cell),
                                 EvaluationFlags::gradients);

      for (const auto q : evaluator_surface.quadrature_point_indices())
        {
          evaluator_surface.submit_normal_derivative(
            evaluator_surface.get_normal_derivative(q), q);
        }

      evaluator_surface.test_and_sum(
        StridedArrayView<Number, n_lanes>(&evaluator.begin_dof_values()[0][v],
                                          dofs_per_cell),
        EvaluationFlags::gradients);

      const unsigned int dofs_per_component = evaluator.dofs_per_component;
      for (unsigned int c = 0; c < components; ++c)
        {
          for (unsigned int i = 0; i < dofs_per_component; ++i)
            deallog << strided_view[c * dofs_per_component + i] << ' ';
          deallog << std::endl;
        }
    }
}
