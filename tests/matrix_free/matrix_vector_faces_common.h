//------------------  matrix_vector_faces_common.h  ------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
//------------------  matrix_vector_faces_common.h  ------------------------


// this is a template for matrix-vector products with face integration (DG
// case, symmetric interior penalty + Nitsche) on different kinds of meshes
// (Cartesian, general, with and without hanging nodes). It also tests the
// multithreading in case it was enabled

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"


// forward declare this function. will be implemented in .cc files
template <int dim, int fe_degree>
void
test();



template <int dim,
          int fe_degree,
          int n_q_points_1d            = fe_degree + 1,
          typename number              = double,
          typename VectorType          = Vector<number>,
          int n_components             = 1,
          typename VectorizedArrayType = VectorizedArray<number>>
class MatrixFreeTest
{
public:
  MatrixFreeTest(const MatrixFree<dim, number, VectorizedArrayType> &data)
    : data(data)
  {}

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = 0;
    data.loop(&MatrixFreeTest::local_apply,
              &MatrixFreeTest::local_apply_face,
              &MatrixFreeTest::local_apply_boundary_face,
              this,
              dst,
              src);
  }

  void
  manual_loop_vmult(VectorType &dst, const VectorType &src) const
  {
    src.update_ghost_values();
    dst = 0;
    local_apply(data, dst, src, std::make_pair(0, data.n_cell_batches()));
    local_apply_face(data,
                     dst,
                     src,
                     std::make_pair(0, data.n_inner_face_batches()));
    local_apply_boundary_face(data,
                              dst,
                              src,
                              std::make_pair(data.n_inner_face_batches(),
                                             data.n_inner_face_batches() +
                                               data.n_boundary_face_batches()));
    local_apply_face(data,
                     dst,
                     src,
                     std::make_pair(data.n_inner_face_batches() +
                                      data.n_boundary_face_batches(),
                                    data.n_inner_face_batches() +
                                      data.n_boundary_face_batches() +
                                      data.n_ghost_inner_face_batches()));
  }

private:
  void
  local_apply(const MatrixFree<dim, number, VectorizedArrayType> &data,
              VectorType                                         &dst,
              const VectorType                                   &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 number,
                 VectorizedArrayType>
      phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);
        phi.evaluate(EvaluationFlags::gradients);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(phi.get_gradient(q), q);
        phi.integrate(EvaluationFlags::gradients);
        phi.distribute_local_to_global(dst);
      }
  }

  void
  local_apply_face(
    const MatrixFree<dim, number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceEvaluation<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     number,
                     VectorizedArrayType>
      fe_eval(data, true);
    FEFaceEvaluation<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     number,
                     VectorizedArrayType>
      fe_eval_neighbor(data, false);
    using value_type =
      typename FEFaceEvaluation<dim,
                                fe_degree,
                                n_q_points_1d,
                                n_components,
                                number,
                                VectorizedArrayType>::value_type;
    const int actual_degree = data.get_dof_handler().get_fe().degree;

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        fe_eval.reinit(face);
        fe_eval_neighbor.reinit(face);

        fe_eval.read_dof_values(src);
        fe_eval.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
        fe_eval_neighbor.read_dof_values(src);
        fe_eval_neighbor.evaluate(EvaluationFlags::values |
                                  EvaluationFlags::gradients);
        VectorizedArrayType sigmaF =
          (std::abs((fe_eval.normal_vector(0) *
                     fe_eval.inverse_jacobian(0))[dim - 1]) +
           std::abs((fe_eval.normal_vector(0) *
                     fe_eval_neighbor.inverse_jacobian(0))[dim - 1])) *
          (number)(std::max(actual_degree, 1) * (actual_degree + 1.0));

        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            value_type average_value =
              (fe_eval.get_value(q) - fe_eval_neighbor.get_value(q)) *
              make_vectorized_array<number, VectorizedArrayType::size()>(0.5);
            value_type average_valgrad =
              fe_eval.get_normal_derivative(q) +
              fe_eval_neighbor.get_normal_derivative(q);
            average_valgrad =
              average_value * sigmaF -
              average_valgrad *
                make_vectorized_array<number, VectorizedArrayType::size()>(0.5);
            fe_eval.submit_normal_derivative(-average_value, q);
            fe_eval_neighbor.submit_normal_derivative(-average_value, q);
            fe_eval.submit_value(average_valgrad, q);
            fe_eval_neighbor.submit_value(-average_valgrad, q);
          }
        fe_eval.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
        fe_eval.distribute_local_to_global(dst);
        fe_eval_neighbor.integrate(EvaluationFlags::values |
                                   EvaluationFlags::gradients);
        fe_eval_neighbor.distribute_local_to_global(dst);
      }
  }


  void
  local_apply_boundary_face(
    const MatrixFree<dim, number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceEvaluation<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     number,
                     VectorizedArrayType>
      fe_eval(data, true);
    using value_type =
      typename FEFaceEvaluation<dim,
                                fe_degree,
                                n_q_points_1d,
                                n_components,
                                number,
                                VectorizedArrayType>::value_type;
    const int actual_degree = data.get_dof_handler().get_fe().degree;
    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        fe_eval.reinit(face);
        fe_eval.read_dof_values(src);
        fe_eval.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
        VectorizedArrayType sigmaF =
          2.0 *
          std::abs(
            (fe_eval.normal_vector(0) * fe_eval.inverse_jacobian(0))[dim - 1]) *
          number(std::max(actual_degree, 1) * (actual_degree + 1.0));

        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            value_type average_value   = fe_eval.get_value(q);
            value_type average_valgrad = -fe_eval.get_normal_derivative(q);
            average_valgrad += average_value * sigmaF;
            fe_eval.submit_normal_derivative(-average_value, q);
            fe_eval.submit_value(average_valgrad, q);
          }

        fe_eval.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
        fe_eval.distribute_local_to_global(dst);
      }
  }

  const MatrixFree<dim, number, VectorizedArrayType> &data;
};



// A variant class used in some of the tests that goes through the combined
// vector-access/evaluate routines that typically give better performance
template <int dim,
          int fe_degree,
          int n_q_points_1d            = fe_degree + 1,
          typename number              = double,
          typename VectorType          = Vector<number>,
          int n_components             = 1,
          typename VectorizedArrayType = VectorizedArray<number>>
class MatrixFreeVariant
{
public:
  MatrixFreeVariant(const MatrixFree<dim, number, VectorizedArrayType> &data,
                    const bool         zero_within_loop       = true,
                    const unsigned int start_vector_component = 0)
    : data(data)
    , zero_within_loop(zero_within_loop)
    , start_vector_component(start_vector_component)
  {}

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    if (!zero_within_loop)
      dst = 0;
    data.loop(&MatrixFreeVariant::local_apply,
              &MatrixFreeVariant::local_apply_face,
              &MatrixFreeVariant::local_apply_boundary_face,
              this,
              dst,
              src,
              zero_within_loop,
              MatrixFree<dim, number, VectorizedArrayType>::DataAccessOnFaces::
                gradients,
              MatrixFree<dim, number, VectorizedArrayType>::DataAccessOnFaces::
                gradients);
  }

  void
  vmult_add(VectorType &dst, const VectorType &src) const
  {
    data.loop(&MatrixFreeVariant::local_apply,
              &MatrixFreeVariant::local_apply_face,
              &MatrixFreeVariant::local_apply_boundary_face,
              this,
              dst,
              src,
              false,
              MatrixFree<dim, number, VectorizedArrayType>::DataAccessOnFaces::
                gradients,
              MatrixFree<dim, number, VectorizedArrayType>::DataAccessOnFaces::
                gradients);
  }

private:
  void
  local_apply(const MatrixFree<dim, number, VectorizedArrayType> &data,
              VectorType                                         &dst,
              const VectorType                                   &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 number,
                 VectorizedArrayType>
      phi(data, 0, 0, start_vector_component);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(src, EvaluationFlags::gradients);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(phi.get_gradient(q), q);
        phi.integrate_scatter(EvaluationFlags::gradients, dst);
      }
  }

  void
  local_apply_face(
    const MatrixFree<dim, number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceEvaluation<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     number,
                     VectorizedArrayType>
      fe_eval(data, true, 0, 0, start_vector_component);
    FEFaceEvaluation<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     number,
                     VectorizedArrayType>
      fe_eval_neighbor(data, false, 0, 0, start_vector_component);
    using value_type =
      typename FEFaceEvaluation<dim,
                                fe_degree,
                                n_q_points_1d,
                                n_components,
                                number,
                                VectorizedArrayType>::value_type;
    const int actual_degree = data.get_dof_handler().get_fe().degree;

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        fe_eval.reinit(face);
        fe_eval_neighbor.reinit(face);

        fe_eval.gather_evaluate(src,
                                EvaluationFlags::values |
                                  EvaluationFlags::gradients);
        fe_eval_neighbor.gather_evaluate(src,
                                         EvaluationFlags::values |
                                           EvaluationFlags::gradients);

        VectorizedArrayType sigmaF =
          (std::abs((fe_eval.normal_vector(0) *
                     fe_eval.inverse_jacobian(0))[dim - 1]) +
           std::abs((fe_eval.normal_vector(0) *
                     fe_eval_neighbor.inverse_jacobian(0))[dim - 1])) *
          (number)(std::max(actual_degree, 1) * (actual_degree + 1.0));

        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            value_type average_value =
              (fe_eval.get_value(q) - fe_eval_neighbor.get_value(q)) *
              make_vectorized_array<number, VectorizedArrayType::size()>(0.5);
            value_type average_valgrad =
              fe_eval.get_normal_derivative(q) +
              fe_eval_neighbor.get_normal_derivative(q);
            average_valgrad =
              average_value * sigmaF -
              average_valgrad *
                make_vectorized_array<number, VectorizedArrayType::size()>(0.5);
            fe_eval.submit_normal_derivative(-average_value, q);
            fe_eval_neighbor.submit_normal_derivative(-average_value, q);
            fe_eval.submit_value(average_valgrad, q);
            fe_eval_neighbor.submit_value(-average_valgrad, q);
          }
        fe_eval.integrate_scatter(EvaluationFlags::values |
                                    EvaluationFlags::gradients,
                                  dst);
        fe_eval_neighbor.integrate_scatter(EvaluationFlags::values |
                                             EvaluationFlags::gradients,
                                           dst);
      }
  }

  void
  local_apply_boundary_face(
    const MatrixFree<dim, number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceEvaluation<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     number,
                     VectorizedArrayType>
      fe_eval(data, true, 0, 0, start_vector_component);
    using value_type =
      typename FEFaceEvaluation<dim,
                                fe_degree,
                                n_q_points_1d,
                                n_components,
                                number,
                                VectorizedArrayType>::value_type;
    const int actual_degree = data.get_dof_handler().get_fe().degree;
    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        fe_eval.reinit(face);
        fe_eval.gather_evaluate(src,
                                EvaluationFlags::values |
                                  EvaluationFlags::gradients);
        VectorizedArrayType sigmaF =
          std::abs(
            (fe_eval.normal_vector(0) * fe_eval.inverse_jacobian(0))[dim - 1]) *
          number(std::max(actual_degree, 1) * (actual_degree + 1.0)) * 2.0;

        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            value_type average_value   = fe_eval.get_value(q);
            value_type average_valgrad = -fe_eval.get_normal_derivative(q);
            average_valgrad += average_value * sigmaF;
            fe_eval.submit_normal_derivative(-average_value, q);
            fe_eval.submit_value(average_valgrad, q);
          }

        fe_eval.integrate_scatter(EvaluationFlags::values |
                                    EvaluationFlags::gradients,
                                  dst);
      }
  }

  const MatrixFree<dim, number, VectorizedArrayType> &data;
  const bool                                          zero_within_loop;
  const unsigned int                                  start_vector_component;
};



template <int dim, int n_components, typename Number>
Tensor<1, n_components, Tensor<1, dim, Number>>
multiply_by_advection(const Tensor<1, dim, Number>          &advection,
                      const Tensor<1, n_components, Number> &values)
{
  Tensor<1, n_components, Tensor<1, dim, Number>> out;
  for (unsigned int c = 0; c < n_components; ++c)
    for (unsigned int d = 0; d < dim; ++d)
      out[c][d] = advection[d] * values[c];
  return out;
}



template <int dim, typename Number>
Tensor<1, dim, Number>
multiply_by_advection(const Tensor<1, dim, Number> &advection,
                      const Number                 &values)
{
  Tensor<1, dim, Number> out;
  for (unsigned int d = 0; d < dim; ++d)
    out[d] = advection[d] * values;
  return out;
}



// An implementation of matrix-free advection operator
template <int dim,
          int fe_degree,
          int n_q_points_1d            = fe_degree + 1,
          typename number              = double,
          typename VectorType          = Vector<number>,
          int n_components             = 1,
          typename VectorizedArrayType = VectorizedArray<number>>
class MatrixFreeAdvection
{
public:
  MatrixFreeAdvection(const MatrixFree<dim, number, VectorizedArrayType> &data,
                      const bool         zero_within_loop       = true,
                      const unsigned int start_vector_component = 0)
    : data(data)
    , zero_within_loop(zero_within_loop)
    , start_vector_component(start_vector_component)
  {
    for (unsigned int d = 0; d < dim; ++d)
      advection[d] = 0.4 + 0.12 * d;
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    if (!zero_within_loop)
      dst = 0;
    data.loop(
      &MatrixFreeAdvection::local_apply,
      &MatrixFreeAdvection::local_apply_face,
      &MatrixFreeAdvection::local_apply_boundary_face,
      this,
      dst,
      src,
      zero_within_loop,
      MatrixFree<dim, number, VectorizedArrayType>::DataAccessOnFaces::values,
      MatrixFree<dim, number, VectorizedArrayType>::DataAccessOnFaces::values);
  }

  void
  vmult_add(VectorType &dst, const VectorType &src) const
  {
    data.loop(
      &MatrixFreeAdvection::local_apply,
      &MatrixFreeAdvection::local_apply_face,
      &MatrixFreeAdvection::local_apply_boundary_face,
      this,
      dst,
      src,
      false,
      MatrixFree<dim, number, VectorizedArrayType>::DataAccessOnFaces::values,
      MatrixFree<dim, number, VectorizedArrayType>::DataAccessOnFaces::values);
  }

private:
  void
  local_apply(const MatrixFree<dim, number, VectorizedArrayType> &data,
              VectorType                                         &dst,
              const VectorType                                   &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 number,
                 VectorizedArrayType>
      phi(data, 0, 0, start_vector_component);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(src, EvaluationFlags::values);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(multiply_by_advection(advection,
                                                    phi.get_value(q)),
                              q);
        phi.integrate_scatter(EvaluationFlags::gradients, dst);
      }
  }

  void
  local_apply_face(
    const MatrixFree<dim, number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceEvaluation<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     number,
                     VectorizedArrayType>
      phi_m(data, true, 0, 0, start_vector_component);
    FEFaceEvaluation<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     number,
                     VectorizedArrayType>
      phi_p(data, false, 0, 0, start_vector_component);
    using value_type =
      typename FEFaceEvaluation<dim,
                                fe_degree,
                                n_q_points_1d,
                                n_components,
                                number,
                                VectorizedArrayType>::value_type;

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi_m.reinit(face);
        phi_m.gather_evaluate(src, EvaluationFlags::values);
        phi_p.reinit(face);
        phi_p.gather_evaluate(src, EvaluationFlags::values);

        for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
          {
            value_type u_minus = phi_m.get_value(q),
                       u_plus  = phi_p.get_value(q);
            const VectorizedArrayType normal_times_advection =
              advection * phi_m.normal_vector(q);
            const value_type flux_times_normal =
              make_vectorized_array<number, VectorizedArrayType::size()>(0.5) *
              ((u_minus + u_plus) * normal_times_advection +
               std::abs(normal_times_advection) * (u_minus - u_plus));
            phi_m.submit_value(-flux_times_normal, q);
            phi_p.submit_value(flux_times_normal, q);
          }

        phi_m.integrate_scatter(EvaluationFlags::values, dst);
        phi_p.integrate_scatter(EvaluationFlags::values, dst);
      }
  }

  void
  local_apply_boundary_face(
    const MatrixFree<dim, number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceEvaluation<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     number,
                     VectorizedArrayType>
      fe_eval(data, true, 0, 0, start_vector_component);
    using value_type =
      typename FEFaceEvaluation<dim,
                                fe_degree,
                                n_q_points_1d,
                                n_components,
                                number,
                                VectorizedArrayType>::value_type;
    value_type u_plus = {};
    for (unsigned int d = 0; d < n_components; ++d)
      u_plus[d] = 1.3;

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        fe_eval.reinit(face);
        fe_eval.gather_evaluate(src, EvaluationFlags::values);

        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            value_type                u_minus = fe_eval.get_value(q);
            const VectorizedArrayType normal_times_advection =
              advection * fe_eval.normal_vector(q);
            const value_type flux_times_normal =
              make_vectorized_array<number, VectorizedArrayType::size()>(0.5) *
              ((u_minus + u_plus) * normal_times_advection +
               std::abs(normal_times_advection) * (u_minus - u_plus));
            fe_eval.submit_value(-flux_times_normal, q);
          }

        fe_eval.integrate_scatter(EvaluationFlags::values, dst);
      }
  }

  const MatrixFree<dim, number, VectorizedArrayType> &data;
  const bool                                          zero_within_loop;
  const unsigned int                                  start_vector_component;
  Tensor<1, dim, VectorizedArrayType>                 advection;
};



// Helper integrators replacing LocalIntegrators::Laplace

template <int dim, typename FEType>
void
assemble_cell_matrix(FullMatrix<double> &M, const FEType &fe)
{
  // assume M has been resized by caller; just zero and assemble
  M                    = 0;
  const unsigned int n = fe.get_fe().dofs_per_cell;
  Assert(M.m() == (int)n && M.n() == (int)n,
         ExcMessage("matrix size mismatch"));
  for (unsigned int q = 0; q < fe.n_quadrature_points; ++q)
    {
      const double w = fe.JxW(q);
      for (unsigned int i = 0; i < n; ++i)
        {
          const Tensor<1, dim> gi = fe.shape_grad(i, q);
          for (unsigned int j = 0; j < n; ++j)
            M(i, j) += (gi * fe.shape_grad(j, q)) * w;
        }
    }
}


template <int dim, typename FEFM, typename FEFP>
void
assemble_interior_face_ip(FullMatrix<double> &Mcc,
                          FullMatrix<double> &Mcn,
                          FullMatrix<double> &Mnc,
                          FullMatrix<double> &Mnn,
                          const FEFM         &fe_m,
                          const FEFP         &fe_p,
                          const double        penalty)
{
  Mcc = 0;
  Mcn = 0;
  Mnc = 0;
  Mnn = 0;

  const unsigned int nm = fe_m.get_fe().dofs_per_cell;
  const unsigned int np = fe_p.get_fe().dofs_per_cell;
  Assert(Mcc.m() == (int)nm && Mcc.n() == (int)nm, ExcMessage("Mcc size"));
  Assert(Mcn.m() == (int)nm && Mcn.n() == (int)np, ExcMessage("Mcn size"));
  Assert(Mnc.m() == (int)np && Mnc.n() == (int)nm, ExcMessage("Mnc size"));
  Assert(Mnn.m() == (int)np && Mnn.n() == (int)np, ExcMessage("Mnn size"));

  // first pass: use minus side normal
  for (unsigned int q = 0; q < fe_m.n_quadrature_points; ++q)
    {
      const double         w    = fe_m.JxW(q);
      const Tensor<1, dim> nvec = fe_m.normal_vector(q);

      for (unsigned int i = 0; i < nm; ++i)
        {
          const double phi_i_m    = fe_m.shape_value(i, q);
          const double dn_phi_i_m = fe_m.shape_grad(i, q) * nvec;
          for (unsigned int j = 0; j < nm; ++j)
            {
              const double phi_j_m    = fe_m.shape_value(j, q);
              const double dn_phi_j_m = fe_m.shape_grad(j, q) * nvec;
              Mcc(i, j) +=
                (penalty * phi_i_m * phi_j_m -
                 0.5 * (dn_phi_j_m * phi_i_m + dn_phi_i_m * phi_j_m)) *
                w;
            }
        }

      for (unsigned int i = 0; i < nm; ++i)
        for (unsigned int j = 0; j < np; ++j)
          {
            const double phi_i_m    = fe_m.shape_value(i, q);
            const double dn_phi_i_m = fe_m.shape_grad(i, q) * nvec;
            const double phi_j_p    = fe_p.shape_value(j, q);
            const double dn_phi_j_p =
              fe_p.shape_grad(j, q) *
              nvec; // note: use same normal direction for consistency

            Mcn(i, j) += (-penalty * phi_i_m * phi_j_p +
                          0.5 * (dn_phi_j_p * phi_i_m - dn_phi_i_m * phi_j_p)) *
                         w;
            Mnc(j, i) += (-penalty * phi_j_p * phi_i_m +
                          0.5 * (dn_phi_i_m * phi_j_p - dn_phi_j_p * phi_i_m)) *
                         w;
          }
    }

  // second pass: fill Mnn correctly using plus side normal
  for (unsigned int q = 0; q < fe_p.n_quadrature_points; ++q)
    {
      const double         w    = fe_p.JxW(q);
      const Tensor<1, dim> nvec = fe_p.normal_vector(q);
      for (unsigned int i = 0; i < np; ++i)
        for (unsigned int j = 0; j < np; ++j)
          {
            const double phi_i_p    = fe_p.shape_value(i, q);
            const double dn_phi_i_p = fe_p.shape_grad(i, q) * nvec;
            const double phi_j_p    = fe_p.shape_value(j, q);
            const double dn_phi_j_p = fe_p.shape_grad(j, q) * nvec;
            Mnn(i, j) += (penalty * phi_i_p * phi_j_p -
                          0.5 * (dn_phi_j_p * phi_i_p + dn_phi_i_p * phi_j_p)) *
                         w;
          }
    }
}


template <int dim, typename FEType>
void
assemble_boundary_nitsche(FullMatrix<double> &Mcc,
                          const FEType       &fe,
                          const double        penalty)
{
  Mcc                  = 0;
  const unsigned int n = fe.get_fe().dofs_per_cell;
  Assert(Mcc.m() == (int)n && Mcc.n() == (int)n, ExcMessage("Mcc size"));

  for (unsigned int q = 0; q < fe.n_quadrature_points; ++q)
    {
      const double         w    = fe.JxW(q);
      const Tensor<1, dim> nvec = fe.normal_vector(q);
      for (unsigned int i = 0; i < n; ++i)
        {
          const double phi_i    = fe.shape_value(i, q);
          const double dn_phi_i = fe.shape_grad(i, q) * nvec;
          for (unsigned int j = 0; j < n; ++j)
            {
              const double phi_j    = fe.shape_value(j, q);
              const double dn_phi_j = fe.shape_grad(j, q) * nvec;
              Mcc(i, j) += (penalty * phi_i * phi_j -
                            0.5 * (dn_phi_j * phi_i + dn_phi_i * phi_j)) *
                           w;
            }
        }
    }
}


// Reference solution created with MeshWorker but implemented via mesh_loop
template <int dim>
class MatrixIntegrator
{
public:
  MatrixIntegrator()
    : matrix(nullptr)
  {}

  MatrixIntegrator(SparseMatrix<double> &matrix)
    : matrix(&matrix)
  {}

  // cell worker
  template <typename CellIteratorType>
  void
  cell(const CellIteratorType                &cell,
       MeshWorker::ScratchData<dim>          &scratch,
       MeshWorker::CopyData<4, 1, 2, double> &copy) const
  {
    const auto        &fe   = scratch.reinit(cell);
    const unsigned int dofs = fe.get_fe().dofs_per_cell;
    copy.reinit(0, dofs, dofs);
    copy.local_dof_indices[0].resize(dofs);
    cell->get_active_or_mg_dof_indices(copy.local_dof_indices[0]);
    // compute cell matrix
    assemble_cell_matrix<dim>(copy.matrices[0], fe);
  }

  // face worker
  template <typename CellIteratorType>
  void
  face(const CellIteratorType                &cell1,
       const unsigned int                     face1,
       const unsigned int                     subface1,
       const CellIteratorType                &cell2,
       const unsigned int                     face2,
       const unsigned int                     subface2,
       MeshWorker::ScratchData<dim>          &scratch,
       MeshWorker::CopyData<4, 1, 2, double> &copy) const
  {
    const auto &fe1 = scratch.reinit(cell1, face1);
    const auto &fe2 = scratch.reinit(cell2, face2);

    const unsigned int dofs1 = fe1.get_fe().dofs_per_cell;
    const unsigned int dofs2 = fe2.get_fe().dofs_per_cell;

    // prepare storage: M_cc, M_cn, M_nc, M_nn
    copy.reinit(0, dofs1, dofs1);
    copy.reinit(1, dofs1, dofs2);
    copy.reinit(2, dofs2, dofs1);
    copy.reinit(3, dofs2, dofs2);

    copy.local_dof_indices[0].resize(dofs1);
    copy.local_dof_indices[1].resize(dofs2);
    cell1->get_active_or_mg_dof_indices(copy.local_dof_indices[0]);
    cell2->get_active_or_mg_dof_indices(copy.local_dof_indices[1]);

    // compute penalty similar to previous implementation
    const unsigned int deg = fe1.get_fe().tensor_degree();
    Tensor<2, dim>     inverse_jacobian =
      transpose(fe1.jacobian(0).covariant_form());
    const double normal_volume_fraction1 = std::abs(
      (inverse_jacobian[GeometryInfo<dim>::unit_normal_direction[face1]] *
       fe1.normal_vector(0)));
    inverse_jacobian = transpose(fe2.jacobian(0).covariant_form());
    const double normal_volume_fraction2 = std::abs(
      (inverse_jacobian[GeometryInfo<dim>::unit_normal_direction[face2]] *
       fe2.normal_vector(0)));
    double penalty = 0.5 * (normal_volume_fraction1 + normal_volume_fraction2) *
                     std::max(1U, deg) * (deg + 1.0);

    // assemble interior-face IP/SIPG contributions
    assemble_interior_face_ip<dim>(copy.matrices[0],
                                   copy.matrices[1],
                                   copy.matrices[2],
                                   copy.matrices[3],
                                   fe1,
                                   fe2,
                                   penalty);
  }

  // boundary worker
  template <typename CellIteratorType>
  void
  boundary(const CellIteratorType                &cell,
           const unsigned int                     face,
           MeshWorker::ScratchData<dim>          &scratch,
           MeshWorker::CopyData<4, 1, 2, double> &copy) const
  {
    const auto        &fe   = scratch.reinit(cell, face);
    const unsigned int dofs = fe.get_fe().dofs_per_cell;
    copy.reinit(0, dofs, dofs);
    copy.local_dof_indices[0].resize(dofs);
    cell->get_active_or_mg_dof_indices(copy.local_dof_indices[0]);

    const unsigned int deg = fe.get_fe().tensor_degree();
    Tensor<2, dim>     inverse_jacobian =
      transpose(fe.jacobian(0).covariant_form());
    const double normal_volume_fraction = std::abs(
      (inverse_jacobian[GeometryInfo<dim>::unit_normal_direction[face]] *
       fe.normal_vector(0)));
    double penalty =
      2.0 * normal_volume_fraction * std::max(1U, deg) * (deg + 1.0);

    assemble_boundary_nitsche<dim>(copy.matrices[0], fe, penalty);
  }

  // copier: assemble copy data into global sparse matrix
  void
  copier(const MeshWorker::CopyData<4, 1, 2, double> &copy) const
  {
    // assemble cell matrix if present
    if (copy.matrices[0].m() > 0 && copy.local_dof_indices[0].size() > 0)
      {
        const auto &rows = copy.local_dof_indices[0];
        for (unsigned int i = 0; i < copy.matrices[0].m(); ++i)
          for (unsigned int j = 0; j < copy.matrices[0].n(); ++j)
            if (matrix)
              matrix->add(rows[i], rows[j], copy.matrices[0](i, j));
      }

    // If ip matrices were filled, assemble them too
    // M_cn (1) : rows = local_dof_indices[0], cols = local_dof_indices[1]
    if (copy.matrices[1].m() > 0)
      {
        const auto &rows = copy.local_dof_indices[0];
        const auto &cols = copy.local_dof_indices[1];
        for (unsigned int i = 0; i < copy.matrices[1].m(); ++i)
          for (unsigned int j = 0; j < copy.matrices[1].n(); ++j)
            if (matrix)
              matrix->add(rows[i], cols[j], copy.matrices[1](i, j));
      }
    // M_nc (2)
    if (copy.matrices[2].m() > 0)
      {
        const auto &rows = copy.local_dof_indices[1];
        const auto &cols = copy.local_dof_indices[0];
        for (unsigned int i = 0; i < copy.matrices[2].m(); ++i)
          for (unsigned int j = 0; j < copy.matrices[2].n(); ++j)
            if (matrix)
              matrix->add(rows[i], cols[j], copy.matrices[2](i, j));
      }
    // M_nn (3)
    if (copy.matrices[3].m() > 0)
      {
        const auto &rows = copy.local_dof_indices[1];
        const auto &cols = copy.local_dof_indices[1];
        for (unsigned int i = 0; i < copy.matrices[3].m(); ++i)
          for (unsigned int j = 0; j < copy.matrices[3].n(); ++j)
            if (matrix)
              matrix->add(rows[i], cols[j], copy.matrices[3](i, j));
      }
  }

private:
  SparseMatrix<double> *matrix;
};



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          typename number,
          typename VectorizedArrayType = VectorizedArray<number>>
void
do_test(const DoFHandler<dim>           &dof,
        const AffineConstraints<double> &constraints,
        const bool                       also_test_parallel = false)
{
  if (std::is_same_v<number, float> == true)
    deallog.push("float");

  deallog << "Testing " << dof.get_fe().get_name();
  deallog << std::endl;
  // std::cout << "Number of cells: " <<
  // dof.get_triangulation().n_active_cells() << std::endl; std::cout << "Number
  // of degrees of freedom: " << dof.n_dofs() << std::endl; std::cout << "Number
  // of constraints: " << constraints.n_constraints() << std::endl;

  MappingQ<dim> mapping(dof.get_fe().degree + 1);

  Vector<number> in(dof.n_dofs()), out(dof.n_dofs());
  Vector<number> out_dist(out);

  // Set random seed for reproducibility
  Testing::srand(42);
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      const double entry = Testing::rand() / (double)RAND_MAX;
      in(i)              = entry;
    }

  // assemble sparse matrix with MeshWorker
  SparsityPattern      sparsity;
  SparseMatrix<double> matrix;
  {
    DynamicSparsityPattern d_sparsity(dof.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof, d_sparsity);
    sparsity.copy_from(d_sparsity);
  }
  matrix.reinit(sparsity);
  using ScratchDataType = MeshWorker::ScratchData<dim>;
  using CopyDataType    = MeshWorker::CopyData<4, 1, 2, double>;

  const QGauss<dim>     quad_cell(n_q_points_1d > 0 ? n_q_points_1d :
                                                  dof.get_fe().degree + 1);
  const QGauss<dim - 1> quad_face(n_q_points_1d > 0 ? n_q_points_1d :
                                                      dof.get_fe().degree + 1);
  const UpdateFlags     cell_update_flags =
    update_values | update_gradients | update_jacobians | update_JxW_values;
  const UpdateFlags face_update_flags = update_values | update_gradients |
                                        update_normal_vectors |
                                        update_JxW_values | update_jacobians;

  ScratchDataType scratch(
    dof.get_fe(), quad_cell, cell_update_flags, quad_face, face_update_flags);
  CopyDataType copy;

  MatrixIntegrator<dim> integrator(matrix);

  using CellIt = decltype(dof.begin_active());

  std::function<void(const CellIt &, ScratchDataType &, CopyDataType &)>
    cell_worker = [&](const CellIt &cell, ScratchDataType &s, CopyDataType &c) {
      integrator.cell(cell, s, c);
    };

  std::function<void(const CopyDataType &)> copier_fn =
    [&](const CopyDataType &c) { integrator.copier(c); };

  std::function<
    void(const CellIt &, const unsigned int, ScratchDataType &, CopyDataType &)>
    boundary_worker =
      [&](const CellIt      &cell,
          const unsigned int face,
          ScratchDataType   &s,
          CopyDataType      &c) { integrator.boundary(cell, face, s, c); };

  std::function<void(const CellIt &,
                     const unsigned int,
                     const unsigned int,
                     const CellIt &,
                     const unsigned int,
                     const unsigned int,
                     ScratchDataType &,
                     CopyDataType &)>
    face_worker = [&](const CellIt      &cell1,
                      const unsigned int face1,
                      const unsigned int subface1,
                      const CellIt      &cell2,
                      const unsigned int face2,
                      const unsigned int subface2,
                      ScratchDataType   &s,
                      CopyDataType      &c) {
      integrator.face(cell1, face1, subface1, cell2, face2, subface2, s, c);
    };

  MeshWorker::mesh_loop(dof.begin_active(),
                        dof.end(),
                        cell_worker,
                        copier_fn,
                        scratch,
                        copy,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_own_interior_faces_once |
                          MeshWorker::assemble_boundary_faces,
                        boundary_worker,
                        face_worker);

  matrix.vmult(out, in);

  // zero constrained dofs
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    if (constraints.is_constrained(i))
      out(i) = 0;

  MatrixFree<dim, number, VectorizedArrayType> mf_data;
  const QGauss<1> quad(n_q_points_1d > 0 ? n_q_points_1d :
                                           dof.get_fe().degree + 1);
  typename MatrixFree<dim, number, VectorizedArrayType>::AdditionalData data;
  data.tasks_parallel_scheme =
    MatrixFree<dim, number, VectorizedArrayType>::AdditionalData::none;
  data.tasks_block_size = 3;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);

  mf_data.reinit(mapping, dof, constraints, quad, data);

  MatrixFreeTest<dim,
                 fe_degree,
                 n_q_points_1d,
                 number,
                 Vector<number>,
                 1,
                 VectorizedArrayType>
    mf(mf_data);
  mf.vmult(out_dist, in);

  out_dist -= out;
  const double diff_norm = out_dist.linfty_norm() / out.linfty_norm();
  deallog << "Norm of difference:          " << diff_norm << std::endl;

  if (also_test_parallel)
    {
      mf_data.clear();
      data.tasks_parallel_scheme =
        MatrixFree<dim, number, VectorizedArrayType>::AdditionalData::
          partition_partition;
      mf_data.reinit(mapping, dof, constraints, quad, data);

      MatrixFreeTest<dim,
                     fe_degree,
                     n_q_points_1d,
                     number,
                     Vector<number>,
                     1,
                     VectorizedArrayType>
        mf(mf_data);
      mf.vmult(out_dist, in);
      out_dist -= out;

      const double diff_norm = out_dist.linfty_norm() / out.linfty_norm();
      deallog << "Norm of difference parallel: " << diff_norm << std::endl;
    }
  deallog << std::endl;

  if (std::is_same_v<number, float> == true)
    deallog.pop();
}



#ifdef DEAL_II_WITH_MPI
int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc,
                                            argv,
                                            testing_max_num_threads());
  mpi_initlog();
#else
int
main(int argc, char **argv)
{
  initlog();
#endif
  deallog << std::setprecision(3);

  {
    deallog.push("2d");
    test<2, 1>();
    test<2, 2>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1>();
    test<3, 2>();
    deallog.pop();
  }
}
