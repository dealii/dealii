// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Check that MatrixFree::get_face_active_fe_index() returns 0 without
// hp-adaptivity if the cells are sorted into categories, as
// MatrixFree::get_cell_active_fe_index() does, such that
// MatrixFreeTools::compute_diagonal() and MatrixFreeTools::compute_matrix()
// with face integrals use the shape information of the only finite element.
// The diagonal and the matrix of an SIPG Laplace operator have to match the
// ones computed without categories.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include "../tests.h"


template <int dim>
void
test()
{
  constexpr int fe_degree   = 1;
  constexpr int n_points    = fe_degree + 1;
  using Number              = double;
  using VectorizedArrayType = VectorizedArray<Number>;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FE_DGQ<dim>(fe_degree));

  AffineConstraints<Number> constraints;
  constraints.close();

  const Number sigma = 10.;

  const auto cell_operation = [&](auto &phi) {
    phi.evaluate(EvaluationFlags::gradients);
    for (unsigned int q = 0; q < phi.n_q_points; ++q)
      phi.submit_gradient(phi.get_gradient(q), q);
    phi.integrate(EvaluationFlags::gradients);
  };

  const auto face_operation = [&](auto &phi_m, auto &phi_p) {
    phi_m.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
      {
        const auto jump = phi_m.get_value(q) - phi_p.get_value(q);
        const auto flux = sigma * jump - 0.5 * (phi_m.get_normal_derivative(q) +
                                                phi_p.get_normal_derivative(q));
        phi_m.submit_normal_derivative(-0.5 * jump, q);
        phi_p.submit_normal_derivative(-0.5 * jump, q);
        phi_m.submit_value(flux, q);
        phi_p.submit_value(-flux, q);
      }
    phi_m.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  };

  const auto boundary_operation = [&](auto &phi) {
    phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    for (unsigned int q = 0; q < phi.n_q_points; ++q)
      {
        const auto value = phi.get_value(q);
        phi.submit_normal_derivative(-value, q);
        phi.submit_value(2. * sigma * value - phi.get_normal_derivative(q), q);
      }
    phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  };

  const auto compute = [&](const bool          use_categories,
                           Vector<Number>     &diagonal,
                           FullMatrix<Number> &matrix) {
    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
    data.mapping_update_flags                = update_gradients;
    data.mapping_update_flags_inner_faces    = update_values | update_gradients;
    data.mapping_update_flags_boundary_faces = update_values | update_gradients;

    if (use_categories)
      {
        data.cell_vectorization_category.resize(tria.n_active_cells());
        for (unsigned int i = 0; i < tria.n_active_cells(); ++i)
          data.cell_vectorization_category[i] = 1 + i % 3;
        data.cell_vectorization_categories_strict = true;
      }

    MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
    matrix_free.reinit(
      MappingQ1<dim>(), dof_handler, constraints, QGauss<1>(n_points), data);

    if (use_categories)
      {
        unsigned int max_index = 0;
        for (unsigned int face = 0;
             face < matrix_free.n_inner_face_batches() +
                      matrix_free.n_boundary_face_batches();
             ++face)
          {
            max_index =
              std::max(max_index,
                       matrix_free.get_face_active_fe_index({face, face + 1},
                                                            true));
            if (face < matrix_free.n_inner_face_batches())
              max_index =
                std::max(max_index,
                         matrix_free.get_face_active_fe_index({face, face + 1},
                                                              false));
          }
        deallog << "Maximum active FE index of the faces: " << max_index
                << std::endl;
      }

    matrix_free.initialize_dof_vector(diagonal);
    MatrixFreeTools::compute_diagonal<dim,
                                      fe_degree,
                                      n_points,
                                      1,
                                      Number,
                                      VectorizedArrayType>(matrix_free,
                                                           diagonal,
                                                           cell_operation,
                                                           face_operation,
                                                           boundary_operation);

    matrix.reinit(dof_handler.n_dofs(), dof_handler.n_dofs());
    MatrixFreeTools::
      compute_matrix<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>(
        matrix_free,
        constraints,
        matrix,
        cell_operation,
        face_operation,
        boundary_operation);
  };

  Vector<Number>     diagonal_reference, diagonal;
  FullMatrix<Number> matrix_reference, matrix;
  compute(false, diagonal_reference, matrix_reference);
  compute(true, diagonal, matrix);

  diagonal.add(-1., diagonal_reference);
  matrix.add(-1., matrix_reference);

  deallog << "Diagonal matches: "
          << (diagonal.l2_norm() < 1e-12 * diagonal_reference.l2_norm())
          << std::endl;
  deallog << "Matrix matches: "
          << (matrix.frobenius_norm() <
              1e-12 * matrix_reference.frobenius_norm())
          << std::endl;
}


int
main()
{
  initlog();

  test<2>();
}
