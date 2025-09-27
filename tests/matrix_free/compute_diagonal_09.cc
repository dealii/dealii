// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Similar to compute_diagonal_08 but with different dof_no and quad_no.

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include <vector>

#include "../tests.h"

using namespace dealii;

template <int dim>
void
test()
{
  const int fe_degree       = 3;
  const int n_points        = fe_degree + 1;
  const int n_components    = 1;
  using Number              = double;
  using VectorizedArrayType = VectorizedArray<Number>;
  using VectorType          = Vector<Number>;

  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria);
  tria.refine_global(0);

  const FE_DGQ<dim> fe_q(fe_degree);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe_q);

  AffineConstraints<Number> constraints;
  constraints.close();

  const unsigned int                         dof_index       = 1;
  const std::vector<const DoFHandler<dim> *> dof_handler_vec = {&dof_handler,
                                                                &dof_handler};
  const std::vector<const AffineConstraints<Number> *> constraints_vec = {
    &constraints, &constraints};

  typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients;
  additional_data.mapping_update_flags_inner_faces =
    update_values | update_gradients;
  additional_data.mapping_update_flags_boundary_faces =
    update_values | update_gradients;

  MappingQ<dim> mapping(1);
  QGauss<1>     quad(fe_degree + 1);

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  matrix_free.reinit(
    mapping, dof_handler_vec, constraints_vec, quad, additional_data);

  const auto cell_operation = [&](auto &phi) {
    phi.evaluate(EvaluationFlags::gradients);
    for (unsigned int q = 0; q < phi.n_q_points; ++q)
      phi.submit_gradient(phi.get_gradient(q), q);
    phi.integrate(EvaluationFlags::gradients);
  };

  const auto face_operation = [&](auto &phi_m, auto &phi_p) {
    phi_m.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    VectorizedArrayType sigmaF =
      (std::abs((phi_m.normal_vector(0) * phi_m.inverse_jacobian(0))[dim - 1]) +
       std::abs(
         (phi_m.normal_vector(0) * phi_p.inverse_jacobian(0))[dim - 1])) *
      (Number)(std::max(fe_degree, 1) * (fe_degree + 1.0));

    for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
      {
        VectorizedArrayType average_value =
          (phi_m.get_value(q) - phi_p.get_value(q)) * 0.5;
        VectorizedArrayType average_valgrad =
          phi_m.get_normal_derivative(q) + phi_p.get_normal_derivative(q);
        average_valgrad = average_value * 2. * sigmaF - average_valgrad * 0.5;
        phi_m.submit_normal_derivative(-average_value, q);
        phi_p.submit_normal_derivative(-average_value, q);
        phi_m.submit_value(average_valgrad, q);
        phi_p.submit_value(-average_valgrad, q);
      }
    phi_m.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  };

  const auto boundary_operation = [&](auto &phi_m) {
    phi_m.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    VectorizedArrayType sigmaF =
      std::abs((phi_m.normal_vector(0) * phi_m.inverse_jacobian(0))[dim - 1]) *
      Number(std::max(fe_degree, 1) * (fe_degree + 1.0)) * 2.0;

    for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
      {
        VectorizedArrayType average_value   = phi_m.get_value(q);
        VectorizedArrayType average_valgrad = -phi_m.get_normal_derivative(q);
        average_valgrad += average_value * sigmaF * 2.0;
        phi_m.submit_normal_derivative(-average_value, q);
        phi_m.submit_value(average_valgrad, q);
      }

    phi_m.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  };

  const auto vmult = [&](auto &dst, const auto &src) {
    matrix_free.template loop<VectorType, VectorType>(
      [&](
        const auto &matrix_free, auto &dst, const auto &src, const auto range) {
        FEEvaluation<dim,
                     fe_degree,
                     n_points,
                     n_components,
                     Number,
                     VectorizedArrayType>
          phi(matrix_free, dof_index);

        for (unsigned int cell = range.first; cell < range.second; ++cell)
          {
            phi.reinit(cell);
            phi.read_dof_values(src);
            cell_operation(phi);
            phi.set_dof_values(dst);
          }
      },
      [&](
        const auto &matrix_free, auto &dst, const auto &src, const auto range) {
        FEFaceEvaluation<dim,
                         fe_degree,
                         n_points,
                         n_components,
                         Number,
                         VectorizedArrayType>
          phi_m(matrix_free, true, dof_index);
        FEFaceEvaluation<dim,
                         fe_degree,
                         n_points,
                         n_components,
                         Number,
                         VectorizedArrayType>
          phi_p(matrix_free, false, dof_index);

        for (unsigned int face = range.first; face < range.second; ++face)
          {
            phi_m.reinit(face);
            phi_p.reinit(face);
            phi_m.read_dof_values(src);
            phi_p.read_dof_values(src);
            face_operation(phi_m, phi_p);
            phi_m.distribute_local_to_global(dst);
            phi_p.distribute_local_to_global(dst);
          }
      },
      [&](const auto &matrix_free,
          auto       &dst,
          const auto &src,
          const auto  face_range) {
        FEFaceEvaluation<dim,
                         fe_degree,
                         n_points,
                         n_components,
                         Number,
                         VectorizedArrayType>
          phi_m(matrix_free, true, dof_index);

        for (unsigned int face = face_range.first; face < face_range.second;
             face++)
          {
            phi_m.reinit(face);
            phi_m.read_dof_values(src);
            boundary_operation(phi_m);
            phi_m.distribute_local_to_global(dst);
          }
      },
      dst,
      src,
      true);
  };

  // Compute matrix and diagonal (brute-force)
  FullMatrix<Number> full_matrix(dof_handler.n_dofs(), dof_handler.n_dofs());
  VectorType         diagonal(dof_handler.n_dofs());

  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    {
      VectorType src(dof_handler.n_dofs());
      VectorType dst(dof_handler.n_dofs());

      src[i] = 1.0;

      vmult(dst, src);

      diagonal[i] = dst[i];
      for (unsigned int j = 0; j < dof_handler.n_dofs(); ++j)
        full_matrix[j][i] = dst[j];
    }


  // Compute matrix via MatrixFreeTools
  FullMatrix<Number> full_matrix_mft(dof_handler.n_dofs(),
                                     dof_handler.n_dofs());

  MatrixFreeTools::compute_matrix<dim,
                                  fe_degree,
                                  n_points,
                                  n_components,
                                  Number,
                                  VectorizedArrayType>(matrix_free,
                                                       constraints,
                                                       full_matrix_mft,
                                                       cell_operation,
                                                       face_operation,
                                                       boundary_operation,
                                                       dof_index);

  for (unsigned int i = 0; i < full_matrix.m(); ++i)
    for (unsigned int j = 0; j < full_matrix.n(); ++j)
      if (std::abs(full_matrix[i][j] - full_matrix_mft[i][j]) > 1e-6)
        {
          full_matrix.print_formatted(
            deallog.get_file_stream(), 3, true, 0, "0.000e+00");
          deallog << std::endl;

          full_matrix_mft.print_formatted(
            deallog.get_file_stream(), 3, true, 0, "0.000e+00");
          deallog << std::endl;
        }

  // Compute diagonal via MatrixFreeTools
  VectorType diagonal_mft;
  matrix_free.initialize_dof_vector(diagonal_mft);

  MatrixFreeTools::compute_diagonal<dim,
                                    fe_degree,
                                    n_points,
                                    n_components,
                                    Number,
                                    VectorizedArrayType>(matrix_free,
                                                         diagonal_mft,
                                                         cell_operation,
                                                         face_operation,
                                                         boundary_operation,
                                                         dof_index);

  for (unsigned int i = 0; i < diagonal.size(); ++i)
    if (std::abs(diagonal[i] - diagonal_mft[i]) > 1e-6)
      {
        diagonal.print(deallog.get_file_stream());
        deallog << std::endl;

        diagonal_mft.print(deallog.get_file_stream());
        deallog << std::endl;
      }

  deallog << "OK!" << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test<2>();
}
