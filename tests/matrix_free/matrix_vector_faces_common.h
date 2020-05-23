//------------------  matrix_vector_faces_common.h  ------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

// To generate a reference solution
#include <deal.II/integrators/laplace.h>

#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/loop.h>

#include <iostream>


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

private:
  void
  local_apply(const MatrixFree<dim, number, VectorizedArrayType> &data,
              VectorType &                                        dst,
              const VectorType &                                  src,
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
        phi.evaluate(false, true, false);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(phi.get_gradient(q), q);
        phi.integrate(false, true);
        phi.distribute_local_to_global(dst);
      }
  }

  void
  local_apply_face(
    const MatrixFree<dim, number, VectorizedArrayType> &data,
    VectorType &                                        dst,
    const VectorType &                                  src,
    const std::pair<unsigned int, unsigned int> &       face_range) const
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
    typedef
      typename FEFaceEvaluation<dim,
                                fe_degree,
                                n_q_points_1d,
                                n_components,
                                number,
                                VectorizedArrayType>::value_type value_type;
    const int actual_degree = data.get_dof_handler().get_fe().degree;

    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        fe_eval.reinit(face);
        fe_eval_neighbor.reinit(face);

        fe_eval.read_dof_values(src);
        fe_eval.evaluate(true, true);
        fe_eval_neighbor.read_dof_values(src);
        fe_eval_neighbor.evaluate(true, true);
        VectorizedArrayType sigmaF =
          (std::abs((fe_eval.get_normal_vector(0) *
                     fe_eval.inverse_jacobian(0))[dim - 1]) +
           std::abs((fe_eval.get_normal_vector(0) *
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
        fe_eval.integrate(true, true);
        fe_eval.distribute_local_to_global(dst);
        fe_eval_neighbor.integrate(true, true);
        fe_eval_neighbor.distribute_local_to_global(dst);
      }
  }


  void
  local_apply_boundary_face(
    const MatrixFree<dim, number, VectorizedArrayType> &data,
    VectorType &                                        dst,
    const VectorType &                                  src,
    const std::pair<unsigned int, unsigned int> &       face_range) const
  {
    FEFaceEvaluation<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     number,
                     VectorizedArrayType>
      fe_eval(data, true);
    typedef
      typename FEFaceEvaluation<dim,
                                fe_degree,
                                n_q_points_1d,
                                n_components,
                                number,
                                VectorizedArrayType>::value_type value_type;
    const int actual_degree = data.get_dof_handler().get_fe().degree;
    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        fe_eval.reinit(face);
        fe_eval.read_dof_values(src);
        fe_eval.evaluate(true, true);
        VectorizedArrayType sigmaF =
          2.0 *
          std::abs((fe_eval.get_normal_vector(0) *
                    fe_eval.inverse_jacobian(0))[dim - 1]) *
          number(std::max(actual_degree, 1) * (actual_degree + 1.0));

        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            value_type average_value   = fe_eval.get_value(q);
            value_type average_valgrad = -fe_eval.get_normal_derivative(q);
            average_valgrad += average_value * sigmaF;
            fe_eval.submit_normal_derivative(-average_value, q);
            fe_eval.submit_value(average_valgrad, q);
          }

        fe_eval.integrate(true, true);
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
              VectorType &                                        dst,
              const VectorType &                                  src,
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
        phi.gather_evaluate(src, false, true);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(phi.get_gradient(q), q);
        phi.integrate_scatter(false, true, dst);
      }
  }

  void
  local_apply_face(
    const MatrixFree<dim, number, VectorizedArrayType> &data,
    VectorType &                                        dst,
    const VectorType &                                  src,
    const std::pair<unsigned int, unsigned int> &       face_range) const
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
    typedef
      typename FEFaceEvaluation<dim,
                                fe_degree,
                                n_q_points_1d,
                                n_components,
                                number,
                                VectorizedArrayType>::value_type value_type;
    const int actual_degree = data.get_dof_handler().get_fe().degree;

    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        fe_eval.reinit(face);
        fe_eval_neighbor.reinit(face);

        fe_eval.gather_evaluate(src, true, true);
        fe_eval_neighbor.gather_evaluate(src, true, true);

        VectorizedArrayType sigmaF =
          (std::abs((fe_eval.get_normal_vector(0) *
                     fe_eval.inverse_jacobian(0))[dim - 1]) +
           std::abs((fe_eval.get_normal_vector(0) *
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
        fe_eval.integrate_scatter(true, true, dst);
        fe_eval_neighbor.integrate_scatter(true, true, dst);
      }
  }

  void
  local_apply_boundary_face(
    const MatrixFree<dim, number, VectorizedArrayType> &data,
    VectorType &                                        dst,
    const VectorType &                                  src,
    const std::pair<unsigned int, unsigned int> &       face_range) const
  {
    FEFaceEvaluation<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     number,
                     VectorizedArrayType>
      fe_eval(data, true, 0, 0, start_vector_component);
    typedef
      typename FEFaceEvaluation<dim,
                                fe_degree,
                                n_q_points_1d,
                                n_components,
                                number,
                                VectorizedArrayType>::value_type value_type;
    const int actual_degree = data.get_dof_handler().get_fe().degree;
    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        fe_eval.reinit(face);
        fe_eval.gather_evaluate(src, true, true);
        VectorizedArrayType sigmaF =
          std::abs((fe_eval.get_normal_vector(0) *
                    fe_eval.inverse_jacobian(0))[dim - 1]) *
          number(std::max(actual_degree, 1) * (actual_degree + 1.0)) * 2.0;

        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            value_type average_value   = fe_eval.get_value(q);
            value_type average_valgrad = -fe_eval.get_normal_derivative(q);
            average_valgrad += average_value * sigmaF;
            fe_eval.submit_normal_derivative(-average_value, q);
            fe_eval.submit_value(average_valgrad, q);
          }

        fe_eval.integrate_scatter(true, true, dst);
      }
  }

  const MatrixFree<dim, number, VectorizedArrayType> &data;
  const bool                                          zero_within_loop;
  const unsigned int                                  start_vector_component;
};



template <int dim, int n_components, typename Number>
Tensor<1, n_components, Tensor<1, dim, Number>>
multiply_by_advection(const Tensor<1, dim, Number> &         advection,
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
                      const Number &                values)
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
              VectorType &                                        dst,
              const VectorType &                                  src,
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
        phi.gather_evaluate(src, true, false);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(multiply_by_advection(advection,
                                                    phi.get_value(q)),
                              q);
        phi.integrate_scatter(false, true, dst);
      }
  }

  void
  local_apply_face(
    const MatrixFree<dim, number, VectorizedArrayType> &data,
    VectorType &                                        dst,
    const VectorType &                                  src,
    const std::pair<unsigned int, unsigned int> &       face_range) const
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
    typedef
      typename FEFaceEvaluation<dim,
                                fe_degree,
                                n_q_points_1d,
                                n_components,
                                number,
                                VectorizedArrayType>::value_type value_type;

    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        phi_m.reinit(face);
        phi_m.gather_evaluate(src, true, false);
        phi_p.reinit(face);
        phi_p.gather_evaluate(src, true, false);

        for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
          {
            value_type u_minus = phi_m.get_value(q),
                       u_plus  = phi_p.get_value(q);
            const VectorizedArrayType normal_times_advection =
              advection * phi_m.get_normal_vector(q);
            const value_type flux_times_normal =
              make_vectorized_array<number, VectorizedArrayType::size()>(0.5) *
              ((u_minus + u_plus) * normal_times_advection +
               std::abs(normal_times_advection) * (u_minus - u_plus));
            phi_m.submit_value(-flux_times_normal, q);
            phi_p.submit_value(flux_times_normal, q);
          }

        phi_m.integrate_scatter(true, false, dst);
        phi_p.integrate_scatter(true, false, dst);
      }
  }

  void
  local_apply_boundary_face(
    const MatrixFree<dim, number, VectorizedArrayType> &data,
    VectorType &                                        dst,
    const VectorType &                                  src,
    const std::pair<unsigned int, unsigned int> &       face_range) const
  {
    FEFaceEvaluation<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     number,
                     VectorizedArrayType>
      fe_eval(data, true, 0, 0, start_vector_component);
    typedef
      typename FEFaceEvaluation<dim,
                                fe_degree,
                                n_q_points_1d,
                                n_components,
                                number,
                                VectorizedArrayType>::value_type value_type;
    value_type                                                   u_plus;
    u_plus = make_vectorized_array<number, VectorizedArrayType::size()>(1.3);

    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        fe_eval.reinit(face);
        fe_eval.gather_evaluate(src, true, false);

        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            value_type                u_minus = fe_eval.get_value(q);
            const VectorizedArrayType normal_times_advection =
              advection * fe_eval.get_normal_vector(q);
            const value_type flux_times_normal =
              make_vectorized_array<number, VectorizedArrayType::size()>(0.5) *
              ((u_minus + u_plus) * normal_times_advection +
               std::abs(normal_times_advection) * (u_minus - u_plus));
            fe_eval.submit_value(-flux_times_normal, q);
          }

        fe_eval.integrate_scatter(true, false, dst);
      }
  }

  const MatrixFree<dim, number, VectorizedArrayType> &data;
  const bool                                          zero_within_loop;
  const unsigned int                                  start_vector_component;
  Tensor<1, dim, VectorizedArrayType>                 advection;
};



// Reference solution created with MeshWorker
template <int dim>
class MatrixIntegrator : public MeshWorker::LocalIntegrator<dim>
{
public:
  void
  cell(MeshWorker::DoFInfo<dim> &                 dinfo,
       typename MeshWorker::IntegrationInfo<dim> &info) const;
  void
  boundary(MeshWorker::DoFInfo<dim> &                 dinfo,
           typename MeshWorker::IntegrationInfo<dim> &info) const;
  void
  face(MeshWorker::DoFInfo<dim> &                 dinfo1,
       MeshWorker::DoFInfo<dim> &                 dinfo2,
       typename MeshWorker::IntegrationInfo<dim> &info1,
       typename MeshWorker::IntegrationInfo<dim> &info2) const;
};



template <int dim>
void
MatrixIntegrator<dim>::cell(
  MeshWorker::DoFInfo<dim> &                 dinfo,
  typename MeshWorker::IntegrationInfo<dim> &info) const
{
  LocalIntegrators::Laplace::cell_matrix(dinfo.matrix(0, false).matrix,
                                         info.fe_values());
}



template <int dim>
void
MatrixIntegrator<dim>::face(
  MeshWorker::DoFInfo<dim> &                 dinfo1,
  MeshWorker::DoFInfo<dim> &                 dinfo2,
  typename MeshWorker::IntegrationInfo<dim> &info1,
  typename MeshWorker::IntegrationInfo<dim> &info2) const
{
  const unsigned int deg = info1.fe_values(0).get_fe().tensor_degree();
  // Manually compute penalty parameter instead of using the function
  // compute_penalty because we do it slightly differently on non-Cartesian
  // meshes.
  Tensor<2, dim> inverse_jacobian =
    transpose(info1.fe_values(0).jacobian(0).covariant_form());
  const double normal_volume_fraction1 =
    std::abs((inverse_jacobian
                [GeometryInfo<dim>::unit_normal_direction[dinfo1.face_number]] *
              info1.fe_values(0).normal_vector(0)));
  inverse_jacobian = transpose(info2.fe_values(0).jacobian(0).covariant_form());
  const double normal_volume_fraction2 =
    std::abs((inverse_jacobian
                [GeometryInfo<dim>::unit_normal_direction[dinfo2.face_number]] *
              info1.fe_values(0).normal_vector(0)));
  double penalty = 0.5 * (normal_volume_fraction1 + normal_volume_fraction2) *
                   std::max(1U, deg) * (deg + 1.0);
  LocalIntegrators::Laplace::ip_matrix(dinfo1.matrix(0, false).matrix,
                                       dinfo1.matrix(0, true).matrix,
                                       dinfo2.matrix(0, true).matrix,
                                       dinfo2.matrix(0, false).matrix,
                                       info1.fe_values(0),
                                       info2.fe_values(0),
                                       penalty);
}



template <int dim>
void
MatrixIntegrator<dim>::boundary(
  MeshWorker::DoFInfo<dim> &                 dinfo,
  typename MeshWorker::IntegrationInfo<dim> &info) const
{
  const unsigned int deg = info.fe_values(0).get_fe().tensor_degree();
  Tensor<2, dim>     inverse_jacobian =
    transpose(info.fe_values(0).jacobian(0).covariant_form());
  const double normal_volume_fraction =
    std::abs((inverse_jacobian
                [GeometryInfo<dim>::unit_normal_direction[dinfo.face_number]] *
              info.fe_values(0).normal_vector(0)));
  double penalty = normal_volume_fraction * std::max(1U, deg) * (deg + 1.0);
  LocalIntegrators::Laplace::nitsche_matrix(dinfo.matrix(0, false).matrix,
                                            info.fe_values(0),
                                            penalty);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          typename number,
          typename VectorizedArrayType = VectorizedArray<number>>
void
do_test(const DoFHandler<dim> &          dof,
        const AffineConstraints<double> &constraints,
        const bool                       also_test_parallel = false)
{
  if (std::is_same<number, float>::value == true)
    deallog.push("float");

  deallog << "Testing " << dof.get_fe().get_name();
  deallog << std::endl;
  // std::cout << "Number of cells: " <<
  // dof.get_triangulation().n_active_cells() << std::endl; std::cout << "Number
  // of degrees of freedom: " << dof.n_dofs() << std::endl; std::cout << "Number
  // of constraints: " << constraints.n_constraints() << std::endl;

  MappingQGeneric<dim> mapping(dof.get_fe().degree + 1);

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
  MeshWorker::IntegrationInfoBox<dim> info_box;
  UpdateFlags                         update_flags =
    update_values | update_gradients | update_jacobians;
  info_box.add_update_flags_all(update_flags);
  info_box.initialize_gauss_quadrature(n_q_points_1d,
                                       n_q_points_1d,
                                       n_q_points_1d);
  info_box.initialize(dof.get_fe(), mapping);

  MeshWorker::DoFInfo<dim> dof_info(dof);

  MeshWorker::Assembler::MatrixSimple<SparseMatrix<double>> assembler;
  assembler.initialize(matrix);

  MatrixIntegrator<dim> integrator;
  MeshWorker::integration_loop<dim, dim>(
    dof.begin_active(), dof.end(), dof_info, info_box, integrator, assembler);

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

  if (std::is_same<number, float>::value == true)
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
