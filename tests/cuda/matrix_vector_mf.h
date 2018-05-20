// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// This is a template for matrix-vector products with the Helmholtz equation
// (zero and first derivatives) on different kinds of meshes (Cartesian,
// general, with and without hanging nodes).

#include "../tests.h"

#include <deal.II/matrix_free/cuda_fe_evaluation.h>
#include <deal.II/matrix_free/cuda_matrix_free.h>

#include <deal.II/lac/cuda_vector.h>

template <int dim, int fe_degree, typename Number, int n_q_points_1d>
class HelmholtzOperatorQuad
{
public:
  __device__ void operator()(
    CUDAWrappers::FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number>*
                       fe_eval,
    const unsigned int q_point) const;
};

template <int dim, int fe_degree, typename Number, int n_q_points_1d>
__device__ void HelmholtzOperatorQuad<dim, fe_degree, Number, n_q_points_1d>::
                operator()(
  CUDAWrappers::FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number>* fe_eval,
  const unsigned int                                                    q) const
{
  fe_eval->submit_value(Number(10) * fe_eval->get_value(q), q);
  fe_eval->submit_gradient(fe_eval->get_gradient(q), q);
}

template <int dim, int fe_degree, typename Number, int n_q_points_1d>
class HelmholtzOperator
{
public:
  __device__ void
  operator()(
    const unsigned int                                          cell,
    const typename CUDAWrappers::MatrixFree<dim, Number>::Data* gpu_data,
    CUDAWrappers::SharedData<dim, Number>*                      shared_data,
    const Number*                                               src,
    Number*                                                     dst) const;

  static const unsigned int n_dofs_1d = fe_degree + 1;
  static const unsigned int n_local_dofs
    = dealii::Utilities::fixed_int_power<fe_degree + 1, dim>::value;
  static const unsigned int n_q_points
    = dealii::Utilities::fixed_int_power<n_q_points_1d, dim>::value;
};

template <int dim, int fe_degree, typename Number, int n_q_points_1d>
__device__ void
HelmholtzOperator<dim, fe_degree, Number, n_q_points_1d>::
operator()(const unsigned int                                          cell,
           const typename CUDAWrappers::MatrixFree<dim, Number>::Data* gpu_data,
           CUDAWrappers::SharedData<dim, Number>* shared_data,
           const Number*                          src,
           Number*                                dst) const
{
  CUDAWrappers::FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> fe_eval(
    cell, gpu_data, shared_data);
  fe_eval.read_dof_values(src);
  fe_eval.evaluate(true, true);
  fe_eval.apply_quad_point_operations(
    HelmholtzOperatorQuad<dim, fe_degree, Number, n_q_points_1d>());
  fe_eval.integrate(true, true);
  fe_eval.distribute_local_to_global(dst);
}

template <int dim,
          int fe_degree,
          typename Number,
          int n_q_points_1d = fe_degree + 1>
class MatrixFreeTest
{
public:
  MatrixFreeTest(const CUDAWrappers::MatrixFree<dim, Number>& data_in)
    : data(data_in){};

  void
  vmult(LinearAlgebra::CUDAWrappers::Vector<Number>&       dst,
        const LinearAlgebra::CUDAWrappers::Vector<Number>& src);

private:
  const CUDAWrappers::MatrixFree<dim, Number>& data;
};

template <int dim, int fe_degree, typename Number, int n_q_points_1d>
void
MatrixFreeTest<dim, fe_degree, Number, n_q_points_1d>::vmult(
  LinearAlgebra::CUDAWrappers::Vector<Number>&       dst,
  const LinearAlgebra::CUDAWrappers::Vector<Number>& src)
{
  dst = static_cast<Number>(0.);
  HelmholtzOperator<dim, fe_degree, Number, n_q_points_1d> helmholtz_operator;
  data.cell_loop(helmholtz_operator, src, dst);
}
