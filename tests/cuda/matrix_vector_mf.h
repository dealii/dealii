// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// This is a template for matrix-vector products with the Helmholtz equation
// (zero and first derivatives) on different kinds of meshes (Cartesian,
// general, with and without hanging nodes).

#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/cuda_fe_evaluation.h>
#include <deal.II/matrix_free/cuda_matrix_free.h>

#include "../tests.h"

template <int dim, int fe_degree, typename Number, int n_q_points_1d>
class HelmholtzOperatorQuad
{
public:
  __device__
  HelmholtzOperatorQuad(Number coef)
    : coef(coef)
  {}

  __device__ void
  operator()(
    CUDAWrappers::FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number>
      *fe_eval) const;

private:
  Number coef;
};



template <int dim, int fe_degree, typename Number, int n_q_points_1d>
__device__ void
HelmholtzOperatorQuad<dim, fe_degree, Number, n_q_points_1d>::operator()(
  CUDAWrappers::FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> *fe_eval)
  const
{
  fe_eval->submit_value(coef * fe_eval->get_value());
  fe_eval->submit_gradient(fe_eval->get_gradient());
}



template <int dim, int fe_degree, typename Number, int n_q_points_1d>
class HelmholtzOperator
{
public:
  HelmholtzOperator(Number *coefficient)
    : coef(coefficient)
  {}

  __device__ void
  operator()(
    const unsigned int                                          cell,
    const typename CUDAWrappers::MatrixFree<dim, Number>::Data *gpu_data,
    CUDAWrappers::SharedData<dim, Number> *                     shared_data,
    const Number *                                              src,
    Number *                                                    dst) const;

  static const unsigned int n_dofs_1d = fe_degree + 1;
  static const unsigned int n_local_dofs =
    dealii::Utilities::pow(fe_degree + 1, dim);
  static const unsigned int n_q_points =
    dealii::Utilities::pow(n_q_points_1d, dim);

  Number *coef;
};



template <int dim, int fe_degree, typename Number, int n_q_points_1d>
__device__ void
HelmholtzOperator<dim, fe_degree, Number, n_q_points_1d>::operator()(
  const unsigned int                                          cell,
  const typename CUDAWrappers::MatrixFree<dim, Number>::Data *gpu_data,
  CUDAWrappers::SharedData<dim, Number> *                     shared_data,
  const Number *                                              src,
  Number *                                                    dst) const
{
  const unsigned int pos = CUDAWrappers::local_q_point_id<dim, Number>(
    cell, gpu_data, n_dofs_1d, n_q_points);

  CUDAWrappers::FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> fe_eval(
    cell, gpu_data, shared_data);
  fe_eval.read_dof_values(src);
  fe_eval.evaluate(true, true);
  fe_eval.apply_for_each_quad_point(
    HelmholtzOperatorQuad<dim, fe_degree, Number, n_q_points_1d>(coef[pos]));
  fe_eval.integrate(true, true);
  fe_eval.distribute_local_to_global(dst);
}



template <int dim, int fe_degree, typename Number, int n_q_points_1d>
class VaryingCoefficientFunctor
{
public:
  VaryingCoefficientFunctor(Number *coefficient)
    : coef(coefficient)
  {}

  __device__ void
  operator()(
    const unsigned int                                          cell,
    const typename CUDAWrappers::MatrixFree<dim, Number>::Data *gpu_data);

  static const unsigned int n_dofs_1d = fe_degree + 1;
  static const unsigned int n_local_dofs =
    dealii::Utilities::pow(fe_degree + 1, dim);
  static const unsigned int n_q_points =
    dealii::Utilities::pow(n_q_points_1d, dim);

private:
  Number *coef;
};



template <int dim, int fe_degree, typename Number, int n_q_points_1d>
__device__ void
VaryingCoefficientFunctor<dim, fe_degree, Number, n_q_points_1d>::operator()(
  const unsigned int                                          cell,
  const typename CUDAWrappers::MatrixFree<dim, Number>::Data *gpu_data)
{
  const unsigned int pos = CUDAWrappers::local_q_point_id<dim, Number>(
    cell, gpu_data, n_dofs_1d, n_q_points);
  const auto q_point =
    CUDAWrappers::get_quadrature_point<dim, Number>(cell, gpu_data, n_dofs_1d);

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
          typename VectorType = LinearAlgebra::CUDAWrappers::Vector<Number>,
          int n_q_points_1d   = fe_degree + 1>
class MatrixFreeTest : public Subscriptor
{
public:
  MatrixFreeTest(const CUDAWrappers::MatrixFree<dim, Number> &data_in,
                 const unsigned int                           size,
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
  const CUDAWrappers::MatrixFree<dim, Number> &data;
  LinearAlgebra::CUDAWrappers::Vector<Number>  coef;
};

template <int dim,
          int fe_degree,
          typename Number,
          typename VectorType,
          int n_q_points_1d>
MatrixFreeTest<dim, fe_degree, Number, VectorType, n_q_points_1d>::
  MatrixFreeTest(const CUDAWrappers::MatrixFree<dim, Number> &data_in,
                 const unsigned int                           size,
                 const bool                                   constant_coeff)
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
  Assert(false, ExcNotImplemented());
  return 0.;
}

template <int dim,
          int fe_degree,
          typename Number,
          typename VectorType,
          int n_q_points_1d>
void
MatrixFreeTest<dim, fe_degree, Number, VectorType, n_q_points_1d>::vmult(
  VectorType &      dst,
  const VectorType &src) const
{
  dst = static_cast<Number>(0.);
  HelmholtzOperator<dim, fe_degree, Number, n_q_points_1d> helmholtz_operator(
    coef.get_values());
  data.cell_loop(helmholtz_operator, src, dst);
  data.copy_constrained_values(src, dst);
}
