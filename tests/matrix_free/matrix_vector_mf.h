// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// this is a template for matrix-vector products with the Helmholtz equation
// (zero and first derivatives) on different kinds of meshes (Cartesian,
// general, with and without hanging nodes). It also tests the multithreading
// in case it was enabled

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"


template <int dim, int fe_degree, typename VectorType, int n_q_points_1d>
void
helmholtz_operator(const MatrixFree<dim, typename VectorType::value_type> &data,
                   VectorType                                             &dst,
                   const VectorType                                       &src,
                   const std::pair<unsigned int, unsigned int> &cell_range)
{
  using Number = typename VectorType::value_type;
  FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> fe_eval(data);
  const unsigned int n_q_points = fe_eval.n_q_points;

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      fe_eval.reinit(cell);
      fe_eval.read_dof_values(src);
      fe_eval.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          fe_eval.submit_value(Number(10) * fe_eval.get_value(q), q);
          fe_eval.submit_gradient(fe_eval.get_gradient(q), q);
        }
      fe_eval.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
      fe_eval.distribute_local_to_global(dst);
    }
}



template <int dim, typename VectorType>
void
helmholtz_operator_no_template(
  const MatrixFree<dim, typename VectorType::value_type> &data,
  VectorType                                             &dst,
  const VectorType                                       &src,
  const std::pair<unsigned int, unsigned int>            &cell_range)
{
  using Number = typename VectorType::value_type;
  FEEvaluation<dim, -1, 0, 1, Number> fe_eval(data, cell_range);
  const unsigned int                  n_q_points = fe_eval.n_q_points;

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      fe_eval.reinit(cell);
      fe_eval.read_dof_values(src);
      fe_eval.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          fe_eval.submit_value(Number(10) * fe_eval.get_value(q), q);
          fe_eval.submit_gradient(fe_eval.get_gradient(q), q);
        }
      fe_eval.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
      fe_eval.distribute_local_to_global(dst);
    }
}



template <int dim, typename VectorType>
void
helmholtz_operator_no_template(
  const MatrixFree<dim, typename VectorType::value_type> &data,
  VectorType                                             &dst,
  const VectorType                                       &src,
  const std::pair<unsigned int, unsigned int>            &cell_range,
  const unsigned int                                      active_fe_index,
  const unsigned int                                      active_quad_index)
{
  using Number = typename VectorType::value_type;
  FEEvaluation<dim, -1, 0, 1, Number> fe_eval(
    data, 0, 0, 0, active_fe_index, active_quad_index);
  const unsigned int n_q_points = fe_eval.n_q_points;

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      fe_eval.reinit(cell);
      fe_eval.read_dof_values(src);
      fe_eval.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          fe_eval.submit_value(Number(10) * fe_eval.get_value(q), q);
          fe_eval.submit_gradient(fe_eval.get_gradient(q), q);
        }
      fe_eval.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
      fe_eval.distribute_local_to_global(dst);
    }
}



template <int dim,
          int fe_degree,
          typename Number,
          typename VectorType = Vector<Number>,
          int n_q_points_1d   = fe_degree + 1>
class MatrixFreeTest
{
public:
  using vector_t = VectorizedArray<Number>;

  MatrixFreeTest(const MatrixFree<dim, Number> &data_in)
    : data(data_in)
  {}

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    if constexpr (std::is_same_v<VectorType, dealii::ArrayView<Number>>)
      for (unsigned int i = 0; i < dst.size(); ++i)
        dst[i] = 0;
    else
      dst = 0;

    const std::function<
      void(const MatrixFree<dim, typename VectorType::value_type> &,
           VectorType &,
           const VectorType &,
           const std::pair<unsigned int, unsigned int> &)>
      wrap = helmholtz_operator<dim, fe_degree, VectorType, n_q_points_1d>;
    data.cell_loop(wrap, dst, src);
  }

private:
  const MatrixFree<dim, Number> &data;
};
