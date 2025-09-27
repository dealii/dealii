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
// on faces produces the same result with FEFaceEvaluation and FEFaceValues.
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
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>

#include "../tests.h"


template <int dim>
void
test_hessians(const unsigned int                             degree,
              const dealii::FE_Poly<dim>                    &fe,
              const dealii::EvaluationFlags::EvaluationFlags evaluation_flags)
{
  using VectorizedArrayType = VectorizedArray<double>;

  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria, Point<dim>(), 1.);
  tria.refine_global(1);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  MappingQGeneric<dim> mapping(2);

  AffineConstraints<double> constraints;
  VectorTools::interpolate_boundary_values(
    mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
  constraints.close();

  // FEFaceEvaluation
  typename MatrixFree<dim, double, VectorizedArrayType>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags_inner_faces =
    update_values | update_gradients | update_hessians;

  MatrixFree<dim, double, VectorizedArrayType> matrix_free;
  matrix_free.reinit(
    mapping, dof_handler, constraints, QGauss<1>(degree + 1), additional_data);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "Working with " << fe.get_name() << " and "
            << dof_handler.n_dofs() << " dofs" << std::endl;

  LinearAlgebra::distributed::Vector<double> src, dst, dst2;
  matrix_free.initialize_dof_vector(src);
  for (auto &v : src)
    v = random_value<double>();

  matrix_free.initialize_dof_vector(dst);
  matrix_free.initialize_dof_vector(dst2);

  // Setup FEFaceValues
  QGauss<dim - 1>                      quad(degree + 1);
  FEFaceValues<dim>                    fe_face_values_m(mapping,
                                     fe,
                                     quad,
                                     update_values | update_gradients |
                                       update_hessians | update_JxW_values);
  FEFaceValues<dim>                    fe_face_values_p(mapping,
                                     fe,
                                     quad,
                                     update_values | update_gradients |
                                       update_hessians | update_JxW_values);
  Vector<double>                       solution_values_m(fe.dofs_per_cell);
  Vector<double>                       solution_values_p(fe.dofs_per_cell);
  std::vector<Tensor<2, dim>>          solution_hessians(quad.size());
  std::vector<Tensor<1, dim>>          solution_gradients(quad.size());
  std::vector<double>                  solution_values(quad.size());
  std::vector<types::global_dof_index> dof_indices_m(fe.dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices_p(fe.dofs_per_cell);
  src.update_ghost_values();
  dst2 = 0;

  matrix_free.template loop<LinearAlgebra::distributed::Vector<double>,
                            LinearAlgebra::distributed::Vector<double>>(
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      (void)matrix_free;
      (void)dst;
      (void)src;
      (void)range;
    },
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      FEFaceEvaluation<dim, -1, 0, 1, double, VectorizedArrayType> fe_eval_m(
        matrix_free, true);
      FEFaceEvaluation<dim, -1, 0, 1, double, VectorizedArrayType> fe_eval_p(
        matrix_free, false);
      for (unsigned int face = range.first; face < range.second; ++face)
        {
          // FEFaceEvaluation
          fe_eval_m.reinit(face);
          fe_eval_m.read_dof_values(src);
          fe_eval_m.evaluate(evaluation_flags);
          for (unsigned int q = 0; q < fe_eval_m.n_q_points; ++q)
            {
              if (evaluation_flags & EvaluationFlags::hessians)
                fe_eval_m.submit_hessian(fe_eval_m.get_hessian(q), q);
              if (evaluation_flags & EvaluationFlags::gradients)
                fe_eval_m.submit_gradient(fe_eval_m.get_gradient(q), q);
              if (evaluation_flags & EvaluationFlags::values)
                fe_eval_m.submit_value(fe_eval_m.get_value(q), q);
            }
          fe_eval_m.integrate(evaluation_flags);
          fe_eval_m.distribute_local_to_global(dst);

          fe_eval_p.reinit(face);
          fe_eval_p.gather_evaluate(src, evaluation_flags);
          for (unsigned int q = 0; q < fe_eval_p.n_q_points; ++q)
            {
              if (evaluation_flags & EvaluationFlags::hessians)
                fe_eval_p.submit_hessian(fe_eval_p.get_hessian(q), q);
              if (evaluation_flags & EvaluationFlags::gradients)
                fe_eval_p.submit_gradient(fe_eval_p.get_gradient(q), q);
              if (evaluation_flags & EvaluationFlags::values)
                fe_eval_p.submit_value(fe_eval_p.get_value(q), q);
            }
          fe_eval_p.integrate(evaluation_flags);
          fe_eval_p.distribute_local_to_global(dst);

          // FEFaceValues
          for (unsigned int v = 0;
               v < matrix_free.n_active_entries_per_face_batch(face);
               ++v)
            {
              const auto face_accessor_inside =
                matrix_free.get_face_iterator(face, v, true);
              const auto face_accessor_outside =
                matrix_free.get_face_iterator(face, v, false);

              fe_face_values_m.reinit(face_accessor_inside.first,
                                      face_accessor_inside.second);

              fe_face_values_p.reinit(face_accessor_outside.first,
                                      face_accessor_outside.second);

              face_accessor_inside.first->get_dof_indices(dof_indices_m);
              face_accessor_outside.first->get_dof_indices(dof_indices_p);
              constraints.get_dof_values(src,
                                         dof_indices_m.begin(),
                                         solution_values_m.begin(),
                                         solution_values_m.end());

              constraints.get_dof_values(src,
                                         dof_indices_p.begin(),
                                         solution_values_p.begin(),
                                         solution_values_p.end());

              for (unsigned int q = 0; q < quad.size(); ++q)
                {
                  double         values = 0.;
                  Tensor<1, dim> gradients;
                  Tensor<2, dim> hessians;
                  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                    {
                      if (evaluation_flags & EvaluationFlags::hessians)
                        hessians += solution_values_m(i) *
                                    fe_face_values_m.shape_hessian(i, q);
                      if (evaluation_flags & EvaluationFlags::gradients)
                        gradients += solution_values_m(i) *
                                     fe_face_values_m.shape_grad(i, q);
                      if (evaluation_flags & EvaluationFlags::values)
                        values += solution_values_m(i) *
                                  fe_face_values_m.shape_value(i, q);
                    }
                  solution_hessians[q]  = hessians * fe_face_values_m.JxW(q);
                  solution_gradients[q] = gradients * fe_face_values_m.JxW(q);
                  solution_values[q]    = values * fe_face_values_m.JxW(q);
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
                          solution_hessians[q],
                          fe_face_values_m.shape_hessian(i, q));
                      if (evaluation_flags & EvaluationFlags::gradients)
                        sum_gradients += solution_gradients[q] *
                                         fe_face_values_m.shape_grad(i, q);
                      if (evaluation_flags & EvaluationFlags::values)
                        sum_values += solution_values[q] *
                                      fe_face_values_m.shape_value(i, q);
                    }
                  solution_values_m(i) =
                    sum_hessians + sum_gradients + sum_values;
                }
              constraints.distribute_local_to_global(solution_values_m,
                                                     dof_indices_m,
                                                     dst2);

              for (unsigned int q = 0; q < quad.size(); ++q)
                {
                  double         values = 0.;
                  Tensor<1, dim> gradients;
                  Tensor<2, dim> hessians;
                  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                    {
                      if (evaluation_flags & EvaluationFlags::hessians)
                        hessians += solution_values_p(i) *
                                    fe_face_values_p.shape_hessian(i, q);
                      if (evaluation_flags & EvaluationFlags::gradients)
                        gradients += solution_values_p(i) *
                                     fe_face_values_p.shape_grad(i, q);
                      if (evaluation_flags & EvaluationFlags::values)
                        values += solution_values_p(i) *
                                  fe_face_values_p.shape_value(i, q);
                    }
                  solution_hessians[q]  = hessians * fe_face_values_p.JxW(q);
                  solution_gradients[q] = gradients * fe_face_values_p.JxW(q);
                  solution_values[q]    = values * fe_face_values_p.JxW(q);
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
                          solution_hessians[q],
                          fe_face_values_p.shape_hessian(i, q));
                      if (evaluation_flags & EvaluationFlags::gradients)
                        sum_gradients += solution_gradients[q] *
                                         fe_face_values_p.shape_grad(i, q);
                      if (evaluation_flags & EvaluationFlags::values)
                        sum_values += solution_values[q] *
                                      fe_face_values_p.shape_value(i, q);
                    }
                  solution_values_p(i) =
                    sum_hessians + sum_gradients + sum_values;
                }
              constraints.distribute_local_to_global(solution_values_p,
                                                     dof_indices_p,
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

int
main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  initlog();
  deallog << std::setprecision(10);

  {
    dealii::EvaluationFlags::EvaluationFlags evaluation_flags =
      EvaluationFlags::hessians;
    deallog << "test_hessians_only" << std::endl;
    for (unsigned int i = 1; i < 3; ++i)
      {
        test_hessians<2>(i, dealii::FE_Q<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_Q<3>(i), evaluation_flags);
        test_hessians<2>(i, dealii::FE_DGQ<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_DGQ<3>(i), evaluation_flags);
      }
  }

  {
    dealii::EvaluationFlags::EvaluationFlags evaluation_flags =
      EvaluationFlags::values | EvaluationFlags::hessians;
    deallog << "test_hessians_with_values" << std::endl;
    for (unsigned int i = 1; i < 3; ++i)
      {
        test_hessians<2>(i, dealii::FE_Q<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_Q<3>(i), evaluation_flags);
        test_hessians<2>(i, dealii::FE_DGQ<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_DGQ<3>(i), evaluation_flags);
      }
  }

  {
    dealii::EvaluationFlags::EvaluationFlags evaluation_flags =
      EvaluationFlags::gradients | EvaluationFlags::hessians;
    deallog << "test_hessians_with_gradients" << std::endl;
    for (unsigned int i = 1; i < 3; ++i)
      {
        test_hessians<2>(i, dealii::FE_Q<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_Q<3>(i), evaluation_flags);
        test_hessians<2>(i, dealii::FE_DGQ<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_DGQ<3>(i), evaluation_flags);
      }
  }

  {
    dealii::EvaluationFlags::EvaluationFlags evaluation_flags =
      EvaluationFlags::values | EvaluationFlags::gradients |
      EvaluationFlags::hessians;
    deallog << "test_hessians_with_gradients_and_values" << std::endl;
    // run the last test also for cubic polynomials
    for (unsigned int i = 1; i < 4; ++i)
      {
        test_hessians<2>(i, dealii::FE_Q<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_Q<3>(i), evaluation_flags);
        test_hessians<2>(i, dealii::FE_DGQ<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_DGQ<3>(i), evaluation_flags);
      }
  }
}
