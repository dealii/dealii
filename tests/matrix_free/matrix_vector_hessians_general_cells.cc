// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test if applying a matrix vector product with a matrix containing hessians
// on faces produces the same result with FEEvaluation and FEValues.
// This is checked for different combinations of EvaluationFlags, FE types and
// polynomial degrees.

#include <deal.II/base/logstream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
test_hessians(const dealii::FE_Poly<dim>                    &fe,
              const dealii::Quadrature<dim>                 &quad,
              const dealii::EvaluationFlags::EvaluationFlags evaluation_flags)
{
  using VectorizedArrayType = VectorizedArray<double>;

  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria);
  tria.refine_global(1);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  MappingQGeneric<dim> mapping(2);

  AffineConstraints<double> constraints;
  VectorTools::interpolate_boundary_values(
    mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
  constraints.close();

  // FEEvaluation
  typename MatrixFree<dim, double, VectorizedArrayType>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags =
    update_values | update_gradients | update_hessians;

  MatrixFree<dim, double, VectorizedArrayType> matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraints, quad, additional_data);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "Working with " << fe.get_name() << " and "
            << dof_handler.n_dofs() << " dofs" << std::endl;

  LinearAlgebra::distributed::Vector<double> src, dst, dst2;
  matrix_free.initialize_dof_vector(src);
  for (auto &v : src)
    v = random_value<double>();

  matrix_free.initialize_dof_vector(dst);
  matrix_free.initialize_dof_vector(dst2);

  // Setup FEValues
  FEValues<dim> fe_values(mapping,
                          fe,
                          quad,
                          update_values | update_gradients | update_hessians |
                            update_JxW_values);

  Vector<double>                       solution_values_local(fe.dofs_per_cell);
  std::vector<Tensor<2, dim>>          solution_hessians(quad.size());
  std::vector<Tensor<1, dim>>          solution_gradients(quad.size());
  std::vector<double>                  solution_values(quad.size());
  std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);
  src.update_ghost_values();
  dst2 = 0;

  matrix_free.template loop<LinearAlgebra::distributed::Vector<double>,
                            LinearAlgebra::distributed::Vector<double>>(
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      FEEvaluation<dim, -1, 0, 1, double, VectorizedArrayType> fe_eval(
        matrix_free);
      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          // FEEvaluation
          fe_eval.reinit(cell);
          fe_eval.gather_evaluate(src, evaluation_flags);
          for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
            {
              if (evaluation_flags & EvaluationFlags::hessians)
                fe_eval.submit_hessian(fe_eval.get_hessian(q), q);
              if (evaluation_flags & EvaluationFlags::gradients)
                fe_eval.submit_gradient(fe_eval.get_gradient(q), q);
              if (evaluation_flags & EvaluationFlags::values)
                fe_eval.submit_value(fe_eval.get_value(q), q);
            }
          fe_eval.integrate(evaluation_flags);
          fe_eval.distribute_local_to_global(dst);

          // FEValues
          for (unsigned int v = 0;
               v < matrix_free.n_active_entries_per_cell_batch(cell);
               ++v)
            {
              const auto cell_iterator = matrix_free.get_cell_iterator(cell, v);

              fe_values.reinit(cell_iterator);

              cell_iterator->get_dof_indices(dof_indices);
              constraints.get_dof_values(src,
                                         dof_indices.begin(),
                                         solution_values_local.begin(),
                                         solution_values_local.end());

              std::vector<Tensor<2, dim>> hessians_function(quad.size());
              fe_values.get_function_hessians(src, hessians_function);

              for (unsigned int q = 0; q < quad.size(); ++q)
                {
                  double         values = 0.;
                  Tensor<1, dim> gradients;
                  Tensor<2, dim> hessians;
                  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                    {
                      if (evaluation_flags & EvaluationFlags::hessians)
                        hessians += solution_values_local(i) *
                                    fe_values.shape_hessian(i, q);
                      if (evaluation_flags & EvaluationFlags::gradients)
                        gradients +=
                          solution_values_local(i) * fe_values.shape_grad(i, q);
                      if (evaluation_flags & EvaluationFlags::values)
                        values += solution_values_local(i) *
                                  fe_values.shape_value(i, q);
                    }
                  solution_hessians[q]  = hessians * fe_values.JxW(q);
                  solution_gradients[q] = gradients * fe_values.JxW(q);
                  solution_values[q]    = values * fe_values.JxW(q);
                }
              for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                {
                  double sum_hessians  = 0.;
                  double sum_gradients = 0.;
                  double sum_values    = 0.;
                  for (unsigned int q = 0; q < quad.size(); ++q)
                    {
                      if (evaluation_flags & EvaluationFlags::hessians)
                        sum_hessians += double_contract<0, 0, 1, 1>(
                          solution_hessians[q], fe_values.shape_hessian(i, q));
                      if (evaluation_flags & EvaluationFlags::gradients)
                        sum_gradients +=
                          solution_gradients[q] * fe_values.shape_grad(i, q);
                      if (evaluation_flags & EvaluationFlags::values)
                        sum_values +=
                          solution_values[q] * fe_values.shape_value(i, q);
                    }
                  solution_values_local(i) =
                    sum_hessians + sum_gradients + sum_values;
                }
              constraints.distribute_local_to_global(solution_values_local,
                                                     dof_indices,
                                                     dst2);
            }
        }
    },
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      (void)matrix_free;
      (void)dst;
      (void)src;
      (void)range;
    },
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      (void)matrix_free;
      (void)dst;
      (void)src;
      (void)range;
    },
    dst,
    src,
    true);

  const unsigned int end_of_print_dst =
    dof_handler.n_dofs() > 9 ? 9 : dof_handler.n_dofs();

  deallog << "dst FEE: ";
  for (unsigned int i = 0; i < end_of_print_dst; ++i)
    deallog << dst[i] << ' ';
  deallog << std::endl;

  deallog << "dst FEV: ";
  for (unsigned int i = 0; i < end_of_print_dst; ++i)
    deallog << dst2[i] << ' ';
  deallog << std::endl;

  // compare solutions of matrix vector product
  {
    dst2 -= dst;

    double error = 0.;
    if (dst.l2_norm() > 0)
      error = dst2.l2_norm() / dst.l2_norm();
    else
      error = dst2.l2_norm();

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      deallog << "FEValues verification: " << error << std::endl << std::endl;
  }
}



void
test_qiterated(dealii::EvaluationFlags::EvaluationFlags evaluation_flags)
{
  deallog << "test_qiterated" << std::endl;
  for (unsigned int i = 1; i < 3; ++i)
    {
      QIterated<2> quad_2(QGauss<1>(i + 1), 2);
      QIterated<3> quad_3(QGauss<1>(i + 1), 2);
      test_hessians<2>(FE_Q<2>(i), quad_2, evaluation_flags);
      if (i < 3)
        test_hessians<3>(FE_Q<3>(i), quad_3, evaluation_flags);
      test_hessians<2>(FE_DGQ<2>(i), quad_2, evaluation_flags);
      if (i < 3)
        test_hessians<3>(FE_DGQ<3>(i), quad_3, evaluation_flags);
    }
}



void
test_qgauss(dealii::EvaluationFlags::EvaluationFlags evaluation_flags)
{
  deallog << "test_qgauss" << std::endl;
  for (unsigned int i = 1; i < 4; ++i)
    {
      QGauss<2> quad_2(i + 1);
      QGauss<3> quad_3(i + 1);
      test_hessians<2>(FE_Q<2>(i), quad_2, evaluation_flags);
      test_hessians<3>(FE_Q<3>(i), quad_3, evaluation_flags);
      test_hessians<2>(FE_DGQ<2>(i), quad_2, evaluation_flags);
      test_hessians<3>(FE_DGQ<3>(i), quad_3, evaluation_flags);
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  initlog();
  deallog << std::setprecision(10);

  {
    EvaluationFlags::EvaluationFlags evaluation_flags =
      EvaluationFlags::hessians;
    deallog << "test_hessians_only" << std::endl;
    test_qgauss(evaluation_flags);
    test_qiterated(evaluation_flags);
  }

  {
    EvaluationFlags::EvaluationFlags evaluation_flags =
      EvaluationFlags::values | EvaluationFlags::hessians;
    deallog << "test_hessians_with_values" << std::endl;
    test_qgauss(evaluation_flags);
    test_qiterated(evaluation_flags);
  }

  {
    EvaluationFlags::EvaluationFlags evaluation_flags =
      EvaluationFlags::gradients | EvaluationFlags::hessians;
    deallog << "test_hessians_with_gradients" << std::endl;
    test_qgauss(evaluation_flags);
    test_qiterated(evaluation_flags);
  }

  {
    EvaluationFlags::EvaluationFlags evaluation_flags =
      EvaluationFlags::values | EvaluationFlags::gradients |
      EvaluationFlags::hessians;
    deallog << "test_hessians_with_gradients_and_values" << std::endl;
    test_qgauss(evaluation_flags);
    test_qiterated(evaluation_flags);
  }
}
