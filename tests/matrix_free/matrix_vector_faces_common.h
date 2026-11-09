// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// this is a template for matrix-vector products with face integration (DG
// case, symmetric interior penalty + Nitsche) on different kinds of meshes
// (Cartesian, general, with and without hanging nodes). It also tests the
// multithreading in case it was enabled

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/meshworker/mesh_loop.h>

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



// Reference solution: assemble the symmetric interior penalty / Nitsche
// matrix with mesh_loop and FEInterfaceValues. It is compared against the
// matrix-free operator below.
template <int dim>
struct MatrixAssemblyScratchData
{
  MatrixAssemblyScratchData(const Mapping<dim>        &mapping,
                            const FiniteElement<dim>  &fe,
                            const Quadrature<dim>     &quadrature,
                            const Quadrature<dim - 1> &face_quadrature)
    : mapping(&mapping)
    , fe(&fe)
    , quadrature(quadrature)
    , face_quadrature(face_quadrature)
    , fe_values(mapping, fe, quadrature, update_gradients | update_JxW_values)
    , fe_face_values(mapping,
                     fe,
                     face_quadrature,
                     update_values | update_gradients | update_JxW_values |
                       update_normal_vectors | update_inverse_jacobians)
    , fe_interface_values(mapping,
                          fe,
                          face_quadrature,
                          update_values | update_gradients | update_JxW_values |
                            update_normal_vectors | update_inverse_jacobians)
  {}

  MatrixAssemblyScratchData(const MatrixAssemblyScratchData &other)
    : MatrixAssemblyScratchData(*other.mapping,
                                *other.fe,
                                other.quadrature,
                                other.face_quadrature)
  {}

  const Mapping<dim>       *mapping;
  const FiniteElement<dim> *fe;
  Quadrature<dim>           quadrature;
  Quadrature<dim - 1>       face_quadrature;
  FEValues<dim>             fe_values;
  FEFaceValues<dim>         fe_face_values;
  FEInterfaceValues<dim>    fe_interface_values;
};



struct MatrixAssemblyCopyData
{
  struct FaceData
  {
    FullMatrix<double>                   matrix;
    std::vector<types::global_dof_index> dof_indices;
  };

  FullMatrix<double>                   cell_matrix;
  std::vector<types::global_dof_index> cell_dof_indices;
  std::vector<FaceData>                face_data;
};



// Penalty parameter of the interior penalty method: the normal part of the
// inverse Jacobian for the face, averaged over the two cells.
template <int dim>
double
normal_volume_fraction(const FEFaceValuesBase<dim> &fe_face_values,
                       const unsigned int           face_no,
                       const Tensor<1, dim>        &normal)
{
  return std::abs(fe_face_values.inverse_jacobian(
                    0)[GeometryInfo<dim>::unit_normal_direction[face_no]] *
                  normal);
}



template <int dim, typename number, typename VectorizedArrayType>
void
assemble_reference_matrix(const DoFHandler<dim> &dof,
                          const Mapping<dim>    &mapping,
                          const unsigned int     n_q_points_1d,
                          SparseMatrix<double>  &matrix)
{
  const FiniteElement<dim> &fe  = dof.get_fe();
  const unsigned int        deg = fe.tensor_degree();
  const unsigned int        n_q = n_q_points_1d > 0 ? n_q_points_1d : deg + 1;
  const unsigned int        n_comp = fe.n_components();

  using Iterator = typename DoFHandler<dim>::active_cell_iterator;

  const auto cell_worker = [&](const Iterator                 &cell,
                               MatrixAssemblyScratchData<dim> &scratch,
                               MatrixAssemblyCopyData         &copy) {
    scratch.fe_values.reinit(cell);
    const unsigned int n_dofs = fe.n_dofs_per_cell();
    copy.cell_matrix.reinit(n_dofs, n_dofs);
    for (unsigned int q = 0; q < scratch.fe_values.n_quadrature_points; ++q)
      for (unsigned int i = 0; i < n_dofs; ++i)
        for (unsigned int j = 0; j < n_dofs; ++j)
          for (unsigned int c = 0; c < n_comp; ++c)
            copy.cell_matrix(i, j) +=
              scratch.fe_values.shape_grad_component(i, q, c) *
              scratch.fe_values.shape_grad_component(j, q, c) *
              scratch.fe_values.JxW(q);
    copy.cell_dof_indices.resize(n_dofs);
    cell->get_active_or_mg_dof_indices(copy.cell_dof_indices);
  };

  // Nitsche boundary terms: 2 gamma u v - d_n u v - u d_n v
  const auto boundary_worker = [&](const Iterator                 &cell,
                                   const unsigned int              face_no,
                                   MatrixAssemblyScratchData<dim> &scratch,
                                   MatrixAssemblyCopyData         &copy) {
    scratch.fe_face_values.reinit(cell, face_no);
    const auto        &fev    = scratch.fe_face_values;
    const unsigned int n_dofs = fe.n_dofs_per_cell();
    const double       penalty =
      normal_volume_fraction<dim>(fev, face_no, fev.normal_vector(0)) *
      std::max(1U, deg) * (deg + 1.0);

    auto &face_data = copy.face_data.emplace_back();
    face_data.matrix.reinit(n_dofs, n_dofs);
    face_data.dof_indices.resize(n_dofs);
    cell->get_active_or_mg_dof_indices(face_data.dof_indices);
    for (unsigned int q = 0; q < fev.n_quadrature_points; ++q)
      for (unsigned int i = 0; i < n_dofs; ++i)
        for (unsigned int j = 0; j < n_dofs; ++j)
          for (unsigned int c = 0; c < n_comp; ++c)
            face_data.matrix(i, j) +=
              (2. * penalty * fev.shape_value_component(i, q, c) *
                 fev.shape_value_component(j, q, c) -
               (fev.normal_vector(q) * fev.shape_grad_component(i, q, c)) *
                 fev.shape_value_component(j, q, c) -
               (fev.normal_vector(q) * fev.shape_grad_component(j, q, c)) *
                 fev.shape_value_component(i, q, c)) *
              fev.JxW(q);
  };

  // Interior penalty: gamma [u][v] - {d_n u}[v] - [u]{d_n v}
  const auto face_worker = [&](const Iterator                 &cell,
                               const unsigned int              f,
                               const unsigned int              sf,
                               const Iterator                 &ncell,
                               const unsigned int              nf,
                               const unsigned int              nsf,
                               MatrixAssemblyScratchData<dim> &scratch,
                               MatrixAssemblyCopyData         &copy) {
    auto &fiv = scratch.fe_interface_values;
    fiv.reinit(cell, f, sf, ncell, nf, nsf);

    const double penalty =
      0.5 *
      (normal_volume_fraction<dim>(fiv.get_fe_face_values(0),
                                   f,
                                   fiv.normal_vector(0)) +
       normal_volume_fraction<dim>(fiv.get_fe_face_values(1),
                                   nf,
                                   fiv.normal_vector(0))) *
      std::max(1U, deg) * (deg + 1.0);

    const unsigned int n_interface_dofs = fiv.n_current_interface_dofs();
    auto              &face_data        = copy.face_data.emplace_back();
    face_data.matrix.reinit(n_interface_dofs, n_interface_dofs);
    face_data.dof_indices = fiv.get_interface_dof_indices();

    const std::vector<Tensor<1, dim>> &normals = fiv.get_normal_vectors();
    for (unsigned int q = 0; q < fiv.n_quadrature_points; ++q)
      for (unsigned int i = 0; i < n_interface_dofs; ++i)
        for (unsigned int j = 0; j < n_interface_dofs; ++j)
          for (unsigned int c = 0; c < n_comp; ++c)
            face_data.matrix(i, j) +=
              (penalty * fiv.jump_in_shape_values(i, q, c) *
                 fiv.jump_in_shape_values(j, q, c) -
               (fiv.average_of_shape_gradients(j, q, c) * normals[q]) *
                 fiv.jump_in_shape_values(i, q, c) -
               (fiv.average_of_shape_gradients(i, q, c) * normals[q]) *
                 fiv.jump_in_shape_values(j, q, c)) *
              fiv.JxW(q);
  };

  const auto copier = [&](const MatrixAssemblyCopyData &copy) {
    for (unsigned int i = 0; i < copy.cell_dof_indices.size(); ++i)
      for (unsigned int j = 0; j < copy.cell_dof_indices.size(); ++j)
        matrix.add(copy.cell_dof_indices[i],
                   copy.cell_dof_indices[j],
                   copy.cell_matrix(i, j));

    for (const auto &face : copy.face_data)
      for (unsigned int i = 0; i < face.dof_indices.size(); ++i)
        for (unsigned int j = 0; j < face.dof_indices.size(); ++j)
          matrix.add(face.dof_indices[i],
                     face.dof_indices[j],
                     face.matrix(i, j));
  };

  MatrixAssemblyScratchData<dim> scratch(mapping,
                                         fe,
                                         QGauss<dim>(n_q),
                                         QGauss<dim - 1>(n_q));
  MatrixAssemblyCopyData         copy;

  MeshWorker::mesh_loop(dof.begin_active(),
                        dof.end(),
                        cell_worker,
                        copier,
                        scratch,
                        copy,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_boundary_faces |
                          MeshWorker::assemble_own_interior_faces_once,
                        boundary_worker,
                        face_worker);
}



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

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

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

  // assemble the reference sparse matrix with mesh_loop
  SparsityPattern      sparsity;
  SparseMatrix<double> matrix;
  {
    DynamicSparsityPattern d_sparsity(dof.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof, d_sparsity);
    sparsity.copy_from(d_sparsity);
  }
  matrix.reinit(sparsity);
  assemble_reference_matrix<dim, number, VectorizedArrayType>(dof,
                                                              mapping,
                                                              n_q_points_1d,
                                                              matrix);

  // compute matrix-vector product
  matrix.vmult(out, in);

  // zero constrained dofs
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    if (constraints.is_constrained(i))
      out(i) = 0;

  // Later on we will still construct a separate MatrixFree instance for
  // performance testing; name it mf_data to match the rest of the test.
  MatrixFree<dim, number, VectorizedArrayType> mf_data;
  const QGauss<1> quad(n_q_points_1d > 0 ? n_q_points_1d :
                                           dof.get_fe().degree + 1);
  typename MatrixFree<dim, number, VectorizedArrayType>::AdditionalData data;
  data.tasks_parallel_scheme =
    MatrixFree<dim, number, VectorizedArrayType>::AdditionalData::none;
  data.tasks_block_size = 3;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values | update_normal_vectors |
     update_contravariant_transformation);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values | update_normal_vectors |
     update_contravariant_transformation);

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
