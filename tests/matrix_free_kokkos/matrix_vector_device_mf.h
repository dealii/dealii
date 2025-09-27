// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// This is a template for matrix-vector products with the Helmholtz equation
// (zero and first derivatives) on different kinds of meshes (Cartesian,
// general, with and without hanging nodes).

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>

#include "../tests.h"

template <int dim, int fe_degree, typename Number, int n_q_points_1d>
class HelmholtzOperatorQuad
{
public:
  DEAL_II_HOST_DEVICE
  HelmholtzOperatorQuad(Number *coef)
    : coef(coef)
  {}

  DEAL_II_HOST_DEVICE void
  operator()(
    Portable::FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> *fe_eval,
    int q_point) const;

  static const unsigned int n_q_points =
    dealii::Utilities::pow(n_q_points_1d, dim);

private:
  const typename Portable::MatrixFree<dim, Number>::Data *gpu_data;
  Number                                                 *coef;
  int                                                     cell;
};



template <int dim, int fe_degree, typename Number, int n_q_points_1d>
DEAL_II_HOST_DEVICE void
HelmholtzOperatorQuad<dim, fe_degree, Number, n_q_points_1d>::operator()(
  Portable::FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> *fe_eval,
  int q_point) const
{
  unsigned int pos = fe_eval->get_matrix_free_data()->local_q_point_id(
    fe_eval->get_current_cell_index(), n_q_points, q_point);
  fe_eval->submit_value(coef[pos] * fe_eval->get_value(q_point), q_point);
  fe_eval->submit_gradient(fe_eval->get_gradient(q_point), q_point);
}



template <int dim, int fe_degree, typename Number, int n_q_points_1d>
class HelmholtzOperator
{
public:
  static const unsigned int n_local_dofs =
    dealii::Utilities::pow(fe_degree + 1, dim);
  static const unsigned int n_q_points =
    dealii::Utilities::pow(n_q_points_1d, dim);

  HelmholtzOperator(Number *coefficient)
    : coef(coefficient)
  {}

  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceVector<Number>                   &src,
             Portable::DeviceVector<Number>                         &dst) const;

  Number *coef;
};



template <int dim, int fe_degree, typename Number, int n_q_points_1d>
DEAL_II_HOST_DEVICE void
HelmholtzOperator<dim, fe_degree, Number, n_q_points_1d>::operator()(
  const typename Portable::MatrixFree<dim, Number>::Data *data,
  const Portable::DeviceVector<Number>                   &src,
  Portable::DeviceVector<Number>                         &dst) const
{
  Portable::FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> fe_eval(
    data);
  fe_eval.read_dof_values(src);
  fe_eval.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
  fe_eval.apply_for_each_quad_point(
    HelmholtzOperatorQuad<dim, fe_degree, Number, n_q_points_1d>(coef));
  fe_eval.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  fe_eval.distribute_local_to_global(dst);
}



template <int dim, int fe_degree, typename Number, int n_q_points_1d>
class VaryingCoefficientFunctor
{
public:
  VaryingCoefficientFunctor(Number *coefficient)
    : coef(coefficient)
  {}

  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *gpu_data,
             const unsigned int                                      cell,
             const unsigned int                                      q) const;

  static const unsigned int n_dofs_1d = fe_degree + 1;
  static const unsigned int n_local_dofs =
    dealii::Utilities::pow(fe_degree + 1, dim);
  static const unsigned int n_q_points =
    dealii::Utilities::pow(n_q_points_1d, dim);

private:
  Number *coef;
};



template <int dim, int fe_degree, typename Number, int n_q_points_1d>
DEAL_II_HOST_DEVICE void
VaryingCoefficientFunctor<dim, fe_degree, Number, n_q_points_1d>::operator()(
  const typename Portable::MatrixFree<dim, Number>::Data *gpu_data,
  const unsigned int                                      cell,
  const unsigned int                                      q) const
{
  const unsigned int pos     = gpu_data->local_q_point_id(cell, n_q_points, q);
  const auto         q_point = gpu_data->get_quadrature_point(cell, q);



  Number p_square = 0.;
  for (unsigned int i = 0; i < dim; ++i)
    {
      Number coord = q_point[i];
      p_square += coord * coord;
    }
  coef[pos] = 10. / (0.05 + 2. * p_square);
}



template <int dim,
          int fe_degree,
          typename Number,
          typename VectorType =
            LinearAlgebra::distributed::Vector<double, MemorySpace::Default>,
          int n_q_points_1d = fe_degree + 1>
class MatrixFreeTest : public EnableObserverPointer
{
public:
  MatrixFreeTest(const Portable::MatrixFree<dim, Number> &data_in,
                 const unsigned int                       size,
                 const bool constant_coeff = true);

  void
  vmult(VectorType &dst, const VectorType &src) const;

  Number
  el(const unsigned int row, const unsigned int col) const;

  types::global_dof_index
  m() const
  {
    return internal_m;
  }

  types::global_dof_index internal_m;

private:
  const Portable::MatrixFree<dim, Number>                         &data;
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> coef;
};

template <int dim,
          int fe_degree,
          typename Number,
          typename VectorType,
          int n_q_points_1d>
MatrixFreeTest<dim, fe_degree, Number, VectorType, n_q_points_1d>::
  MatrixFreeTest(const Portable::MatrixFree<dim, Number> &data_in,
                 const unsigned int                       size,
                 const bool                               constant_coeff)
  : data(data_in)
{
  coef.reinit(size);
  if (constant_coeff)
    {
      coef.add(10.);
    }
  else
    {
      VaryingCoefficientFunctor<dim, fe_degree, Number, n_q_points_1d> functor(
        coef.get_values());
      data.evaluate_coefficients(functor);
    }
}

template <int dim,
          int fe_degree,
          typename Number,
          typename VectorType,
          int n_q_points_1d>
Number
MatrixFreeTest<dim, fe_degree, Number, VectorType, n_q_points_1d>::el(
  const unsigned int row,
  const unsigned int col) const
{
  (void)col;
  DEAL_II_NOT_IMPLEMENTED();
  return 0.;
}

template <int dim,
          int fe_degree,
          typename Number,
          typename VectorType,
          int n_q_points_1d>
void
MatrixFreeTest<dim, fe_degree, Number, VectorType, n_q_points_1d>::vmult(
  VectorType       &dst,
  const VectorType &src) const
{
  dst = static_cast<Number>(0.);
  HelmholtzOperator<dim, fe_degree, Number, n_q_points_1d> helmholtz_operator(
    coef.get_values());
  data.cell_loop(helmholtz_operator, src, dst);
  data.copy_constrained_values(src, dst);
}
