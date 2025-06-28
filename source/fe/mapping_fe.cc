// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/array_view.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_poly.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_internal.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/full_matrix.h>

#include <boost/container/small_vector.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numeric>


DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
MappingFE<dim, spacedim>::InternalData::InternalData(
  const FiniteElement<dim, spacedim> &fe)
  : fe(fe)
  , polynomial_degree(fe.tensor_degree())
  , n_shape_functions(fe.n_dofs_per_cell())
{}



template <int dim, int spacedim>
std::size_t
MappingFE<dim, spacedim>::InternalData::memory_consumption() const
{
  return (
    Mapping<dim, spacedim>::InternalDataBase::memory_consumption() +
    MemoryConsumption::memory_consumption(shape_values) +
    MemoryConsumption::memory_consumption(shape_derivatives) +
    MemoryConsumption::memory_consumption(covariant) +
    MemoryConsumption::memory_consumption(contravariant) +
    MemoryConsumption::memory_consumption(unit_tangentials) +
    MemoryConsumption::memory_consumption(aux) +
    MemoryConsumption::memory_consumption(mapping_support_points) +
    MemoryConsumption::memory_consumption(cell_of_current_support_points) +
    MemoryConsumption::memory_consumption(volume_elements) +
    MemoryConsumption::memory_consumption(polynomial_degree) +
    MemoryConsumption::memory_consumption(n_shape_functions));
}



template <int dim, int spacedim>
void
MappingFE<dim, spacedim>::InternalData::reinit(const UpdateFlags update_flags,
                                               const Quadrature<dim> &q)
{
  // store the flags in the internal data object so we can access them
  // in fill_fe_*_values()
  this->update_each = update_flags;

  const unsigned int n_q_points = q.size();

  if (this->update_each & update_covariant_transformation)
    covariant.resize(n_q_points);

  if (this->update_each & update_contravariant_transformation)
    contravariant.resize(n_q_points);

  if (this->update_each & update_volume_elements)
    volume_elements.resize(n_q_points);

  // see if we need the (transformation) shape function values
  // and/or gradients and resize the necessary arrays
  if (this->update_each & update_quadrature_points)
    shape_values.resize(n_shape_functions * n_q_points);

  if (this->update_each &
      (update_covariant_transformation | update_contravariant_transformation |
       update_JxW_values | update_boundary_forms | update_normal_vectors |
       update_jacobians | update_jacobian_grads | update_inverse_jacobians |
       update_jacobian_pushed_forward_grads | update_jacobian_2nd_derivatives |
       update_jacobian_pushed_forward_2nd_derivatives |
       update_jacobian_3rd_derivatives |
       update_jacobian_pushed_forward_3rd_derivatives))
    shape_derivatives.resize(n_shape_functions * n_q_points);

  if (this->update_each &
      (update_jacobian_grads | update_jacobian_pushed_forward_grads))
    shape_second_derivatives.resize(n_shape_functions * n_q_points);

  if (this->update_each & (update_jacobian_2nd_derivatives |
                           update_jacobian_pushed_forward_2nd_derivatives))
    shape_third_derivatives.resize(n_shape_functions * n_q_points);

  if (this->update_each & (update_jacobian_3rd_derivatives |
                           update_jacobian_pushed_forward_3rd_derivatives))
    shape_fourth_derivatives.resize(n_shape_functions * n_q_points);

  // now also fill the various fields with their correct values
  compute_shape_function_values(q.get_points());

  // copy (projected) quadrature weights
  quadrature_weights = q.get_weights();
}



template <int dim, int spacedim>
void
MappingFE<dim, spacedim>::InternalData::initialize_face(
  const UpdateFlags      update_flags,
  const Quadrature<dim> &q,
  const unsigned int     n_original_q_points)
{
  reinit(update_flags, q);

  if (this->update_each &
      (update_boundary_forms | update_normal_vectors | update_jacobians |
       update_JxW_values | update_inverse_jacobians))
    {
      aux.resize(dim - 1,
                 std::vector<Tensor<1, spacedim>>(n_original_q_points));

      // Compute tangentials to the unit cell.
      const auto reference_cell = this->fe.reference_cell();
      const auto n_faces        = reference_cell.n_faces();

      for (unsigned int i = 0; i < n_faces; ++i)
        {
          unit_tangentials[i].resize(n_original_q_points);
          std::fill(unit_tangentials[i].begin(),
                    unit_tangentials[i].end(),
                    reference_cell.template face_tangent_vector<dim>(i, 0));
          if (dim > 2)
            {
              unit_tangentials[n_faces + i].resize(n_original_q_points);
              std::fill(unit_tangentials[n_faces + i].begin(),
                        unit_tangentials[n_faces + i].end(),
                        reference_cell.template face_tangent_vector<dim>(i, 1));
            }
        }
    }
}



template <int dim, int spacedim>
void
MappingFE<dim, spacedim>::InternalData::compute_shape_function_values(
  const std::vector<Point<dim>> &unit_points)
{
  const auto fe_poly = dynamic_cast<const FE_Poly<dim, spacedim> *>(&this->fe);

  Assert(fe_poly != nullptr, ExcNotImplemented());

  const auto &tensor_pols = fe_poly->get_poly_space();

  const unsigned int n_shape_functions = fe.n_dofs_per_cell();
  const unsigned int n_points          = unit_points.size();

  std::vector<double>         values;
  std::vector<Tensor<1, dim>> grads;
  if (shape_values.size() != 0)
    {
      Assert(shape_values.size() == n_shape_functions * n_points,
             ExcInternalError());
      values.resize(n_shape_functions);
    }
  if (shape_derivatives.size() != 0)
    {
      Assert(shape_derivatives.size() == n_shape_functions * n_points,
             ExcInternalError());
      grads.resize(n_shape_functions);
    }

  std::vector<Tensor<2, dim>> grad2;
  if (shape_second_derivatives.size() != 0)
    {
      Assert(shape_second_derivatives.size() == n_shape_functions * n_points,
             ExcInternalError());
      grad2.resize(n_shape_functions);
    }

  std::vector<Tensor<3, dim>> grad3;
  if (shape_third_derivatives.size() != 0)
    {
      Assert(shape_third_derivatives.size() == n_shape_functions * n_points,
             ExcInternalError());
      grad3.resize(n_shape_functions);
    }

  std::vector<Tensor<4, dim>> grad4;
  if (shape_fourth_derivatives.size() != 0)
    {
      Assert(shape_fourth_derivatives.size() == n_shape_functions * n_points,
             ExcInternalError());
      grad4.resize(n_shape_functions);
    }


  if (shape_values.size() != 0 || shape_derivatives.size() != 0 ||
      shape_second_derivatives.size() != 0 ||
      shape_third_derivatives.size() != 0 ||
      shape_fourth_derivatives.size() != 0)
    for (unsigned int point = 0; point < n_points; ++point)
      {
        tensor_pols.evaluate(
          unit_points[point], values, grads, grad2, grad3, grad4);

        if (shape_values.size() != 0)
          for (unsigned int i = 0; i < n_shape_functions; ++i)
            shape(point, i) = values[i];

        if (shape_derivatives.size() != 0)
          for (unsigned int i = 0; i < n_shape_functions; ++i)
            derivative(point, i) = grads[i];

        if (shape_second_derivatives.size() != 0)
          for (unsigned int i = 0; i < n_shape_functions; ++i)
            second_derivative(point, i) = grad2[i];

        if (shape_third_derivatives.size() != 0)
          for (unsigned int i = 0; i < n_shape_functions; ++i)
            third_derivative(point, i) = grad3[i];

        if (shape_fourth_derivatives.size() != 0)
          for (unsigned int i = 0; i < n_shape_functions; ++i)
            fourth_derivative(point, i) = grad4[i];
      }
}


namespace internal
{
  namespace MappingFEImplementation
  {
    namespace
    {
      /**
       * Compute the locations of quadrature points on the object described by
       * the first argument (and the cell for which the mapping support points
       * have already been set), but only if the update_flags of the @p data
       * argument indicate so.
       */
      template <int dim, int spacedim>
      void
      maybe_compute_q_points(
        const typename QProjector<dim>::DataSetDescriptor              data_set,
        const typename dealii::MappingFE<dim, spacedim>::InternalData &data,
        std::vector<Point<spacedim>> &quadrature_points,
        const unsigned int            n_q_points)
      {
        const UpdateFlags update_flags = data.update_each;

        if (update_flags & update_quadrature_points)
          for (unsigned int point = 0; point < n_q_points; ++point)
            {
              const double   *shape = &data.shape(point + data_set, 0);
              Point<spacedim> result =
                (shape[0] * data.mapping_support_points[0]);
              for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                for (unsigned int i = 0; i < spacedim; ++i)
                  result[i] += shape[k] * data.mapping_support_points[k][i];
              quadrature_points[point] = result;
            }
      }



      /**
       * Update the co- and contravariant matrices as well as their determinant,
       * for the cell
       * described stored in the data object, but only if the update_flags of the @p data
       * argument indicate so.
       *
       * Skip the computation if possible as indicated by the first argument.
       */
      template <int dim, int spacedim>
      void
      maybe_update_Jacobians(
        const CellSimilarity::Similarity cell_similarity,
        const typename dealii::QProjector<dim>::DataSetDescriptor      data_set,
        const typename dealii::MappingFE<dim, spacedim>::InternalData &data,
        const unsigned int n_q_points)
      {
        const UpdateFlags update_flags = data.update_each;

        if (update_flags & update_contravariant_transformation)
          // if the current cell is just a
          // translation of the previous one, no
          // need to recompute jacobians...
          if (cell_similarity != CellSimilarity::translation)
            {
              std::fill(data.contravariant.begin(),
                        data.contravariant.end(),
                        DerivativeForm<1, dim, spacedim>());

              Assert(data.n_shape_functions > 0, ExcInternalError());

              for (unsigned int point = 0; point < n_q_points; ++point)
                {
                  double result[spacedim][dim];

                  // peel away part of sum to avoid zeroing the
                  // entries and adding for the first time
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim; ++j)
                      result[i][j] = data.derivative(point + data_set, 0)[j] *
                                     data.mapping_support_points[0][i];
                  for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim; ++j)
                        result[i][j] +=
                          data.derivative(point + data_set, k)[j] *
                          data.mapping_support_points[k][i];

                  // write result into contravariant data. for
                  // j=dim in the case dim<spacedim, there will
                  // never be any nonzero data that arrives in
                  // here, so it is ok anyway because it was
                  // initialized to zero at the initialization
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim; ++j)
                      data.contravariant[point][i][j] = result[i][j];
                }
            }

        if (update_flags & update_covariant_transformation)
          if (cell_similarity != CellSimilarity::translation)
            {
              for (unsigned int point = 0; point < n_q_points; ++point)
                {
                  data.covariant[point] =
                    (data.contravariant[point]).covariant_form();
                }
            }

        if (update_flags & update_volume_elements)
          if (cell_similarity != CellSimilarity::translation)
            {
              for (unsigned int point = 0; point < n_q_points; ++point)
                data.volume_elements[point] =
                  data.contravariant[point].determinant();
            }
      }

      /**
       * Update the Hessian of the transformation from unit to real cell, the
       * Jacobian gradients.
       *
       * Skip the computation if possible as indicated by the first argument.
       */
      template <int dim, int spacedim>
      void
      maybe_update_jacobian_grads(
        const CellSimilarity::Similarity                  cell_similarity,
        const typename QProjector<dim>::DataSetDescriptor data_set,
        const typename dealii::MappingFE<dim, spacedim>::InternalData &data,
        std::vector<DerivativeForm<2, dim, spacedim>> &jacobian_grads,
        const unsigned int                             n_q_points)
      {
        const UpdateFlags update_flags = data.update_each;
        if (update_flags & update_jacobian_grads)
          {
            AssertIndexRange(n_q_points, jacobian_grads.size() + 1);

            if (cell_similarity != CellSimilarity::translation)
              for (unsigned int point = 0; point < n_q_points; ++point)
                {
                  const Tensor<2, dim> *second =
                    &data.second_derivative(point + data_set, 0);
                  double result[spacedim][dim][dim];
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim; ++j)
                      for (unsigned int l = 0; l < dim; ++l)
                        result[i][j][l] =
                          (second[0][j][l] * data.mapping_support_points[0][i]);
                  for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          result[i][j][l] +=
                            (second[k][j][l] *
                             data.mapping_support_points[k][i]);

                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim; ++j)
                      for (unsigned int l = 0; l < dim; ++l)
                        jacobian_grads[point][i][j][l] = result[i][j][l];
                }
          }
      }

      /**
       * Update the Hessian of the transformation from unit to real cell, the
       * Jacobian gradients, pushed forward to the real cell coordinates.
       *
       * Skip the computation if possible as indicated by the first argument.
       */
      template <int dim, int spacedim>
      void
      maybe_update_jacobian_pushed_forward_grads(
        const CellSimilarity::Similarity                  cell_similarity,
        const typename QProjector<dim>::DataSetDescriptor data_set,
        const typename dealii::MappingFE<dim, spacedim>::InternalData &data,
        std::vector<Tensor<3, spacedim>> &jacobian_pushed_forward_grads,
        const unsigned int                n_q_points)
      {
        const UpdateFlags update_flags = data.update_each;
        if (update_flags & update_jacobian_pushed_forward_grads)
          {
            AssertIndexRange(n_q_points,
                             jacobian_pushed_forward_grads.size() + 1);

            if (cell_similarity != CellSimilarity::translation)
              {
                double tmp[spacedim][spacedim][spacedim];
                for (unsigned int point = 0; point < n_q_points; ++point)
                  {
                    const Tensor<2, dim> *second =
                      &data.second_derivative(point + data_set, 0);
                    double result[spacedim][dim][dim];
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          result[i][j][l] = (second[0][j][l] *
                                             data.mapping_support_points[0][i]);
                    for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                      for (unsigned int i = 0; i < spacedim; ++i)
                        for (unsigned int j = 0; j < dim; ++j)
                          for (unsigned int l = 0; l < dim; ++l)
                            result[i][j][l] +=
                              (second[k][j][l] *
                               data.mapping_support_points[k][i]);

                    // first push forward the j-components
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < spacedim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          {
                            tmp[i][j][l] =
                              result[i][0][l] * data.covariant[point][j][0];
                            for (unsigned int jr = 1; jr < dim; ++jr)
                              {
                                tmp[i][j][l] += result[i][jr][l] *
                                                data.covariant[point][j][jr];
                              }
                          }

                    // now, pushing forward the l-components
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < spacedim; ++j)
                        for (unsigned int l = 0; l < spacedim; ++l)
                          {
                            jacobian_pushed_forward_grads[point][i][j][l] =
                              tmp[i][j][0] * data.covariant[point][l][0];
                            for (unsigned int lr = 1; lr < dim; ++lr)
                              {
                                jacobian_pushed_forward_grads[point][i][j][l] +=
                                  tmp[i][j][lr] * data.covariant[point][l][lr];
                              }
                          }
                  }
              }
          }
      }

      /**
       * Update the third derivatives of the transformation from unit to real
       * cell, the Jacobian hessians.
       *
       * Skip the computation if possible as indicated by the first argument.
       */
      template <int dim, int spacedim>
      void
      maybe_update_jacobian_2nd_derivatives(
        const CellSimilarity::Similarity                  cell_similarity,
        const typename QProjector<dim>::DataSetDescriptor data_set,
        const typename dealii::MappingFE<dim, spacedim>::InternalData &data,
        std::vector<DerivativeForm<3, dim, spacedim>> &jacobian_2nd_derivatives,
        const unsigned int                             n_q_points)
      {
        const UpdateFlags update_flags = data.update_each;
        if (update_flags & update_jacobian_2nd_derivatives)
          {
            AssertIndexRange(n_q_points, jacobian_2nd_derivatives.size() + 1);

            if (cell_similarity != CellSimilarity::translation)
              {
                for (unsigned int point = 0; point < n_q_points; ++point)
                  {
                    const Tensor<3, dim> *third =
                      &data.third_derivative(point + data_set, 0);
                    double result[spacedim][dim][dim][dim];
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          for (unsigned int m = 0; m < dim; ++m)
                            result[i][j][l][m] =
                              (third[0][j][l][m] *
                               data.mapping_support_points[0][i]);
                    for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                      for (unsigned int i = 0; i < spacedim; ++i)
                        for (unsigned int j = 0; j < dim; ++j)
                          for (unsigned int l = 0; l < dim; ++l)
                            for (unsigned int m = 0; m < dim; ++m)
                              result[i][j][l][m] +=
                                (third[k][j][l][m] *
                                 data.mapping_support_points[k][i]);

                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          for (unsigned int m = 0; m < dim; ++m)
                            jacobian_2nd_derivatives[point][i][j][l][m] =
                              result[i][j][l][m];
                  }
              }
          }
      }

      /**
       * Update the Hessian of the Hessian of the transformation from unit
       * to real cell, the Jacobian Hessian gradients, pushed forward to the
       * real cell coordinates.
       *
       * Skip the computation if possible as indicated by the first argument.
       */
      template <int dim, int spacedim>
      void
      maybe_update_jacobian_pushed_forward_2nd_derivatives(
        const CellSimilarity::Similarity                  cell_similarity,
        const typename QProjector<dim>::DataSetDescriptor data_set,
        const typename dealii::MappingFE<dim, spacedim>::InternalData &data,
        std::vector<Tensor<4, spacedim>>
                          &jacobian_pushed_forward_2nd_derivatives,
        const unsigned int n_q_points)
      {
        const UpdateFlags update_flags = data.update_each;
        if (update_flags & update_jacobian_pushed_forward_2nd_derivatives)
          {
            AssertIndexRange(n_q_points,
                             jacobian_pushed_forward_2nd_derivatives.size() +
                               1);

            if (cell_similarity != CellSimilarity::translation)
              {
                double tmp[spacedim][spacedim][spacedim][spacedim];
                for (unsigned int point = 0; point < n_q_points; ++point)
                  {
                    const Tensor<3, dim> *third =
                      &data.third_derivative(point + data_set, 0);
                    double result[spacedim][dim][dim][dim];
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          for (unsigned int m = 0; m < dim; ++m)
                            result[i][j][l][m] =
                              (third[0][j][l][m] *
                               data.mapping_support_points[0][i]);
                    for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                      for (unsigned int i = 0; i < spacedim; ++i)
                        for (unsigned int j = 0; j < dim; ++j)
                          for (unsigned int l = 0; l < dim; ++l)
                            for (unsigned int m = 0; m < dim; ++m)
                              result[i][j][l][m] +=
                                (third[k][j][l][m] *
                                 data.mapping_support_points[k][i]);

                    // push forward the j-coordinate
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < spacedim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          for (unsigned int m = 0; m < dim; ++m)
                            {
                              jacobian_pushed_forward_2nd_derivatives
                                [point][i][j][l][m] =
                                  result[i][0][l][m] *
                                  data.covariant[point][j][0];
                              for (unsigned int jr = 1; jr < dim; ++jr)
                                jacobian_pushed_forward_2nd_derivatives[point]
                                                                       [i][j][l]
                                                                       [m] +=
                                  result[i][jr][l][m] *
                                  data.covariant[point][j][jr];
                            }

                    // push forward the l-coordinate
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < spacedim; ++j)
                        for (unsigned int l = 0; l < spacedim; ++l)
                          for (unsigned int m = 0; m < dim; ++m)
                            {
                              tmp[i][j][l][m] =
                                jacobian_pushed_forward_2nd_derivatives[point]
                                                                       [i][j][0]
                                                                       [m] *
                                data.covariant[point][l][0];
                              for (unsigned int lr = 1; lr < dim; ++lr)
                                tmp[i][j][l][m] +=
                                  jacobian_pushed_forward_2nd_derivatives
                                    [point][i][j][lr][m] *
                                  data.covariant[point][l][lr];
                            }

                    // push forward the m-coordinate
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < spacedim; ++j)
                        for (unsigned int l = 0; l < spacedim; ++l)
                          for (unsigned int m = 0; m < spacedim; ++m)
                            {
                              jacobian_pushed_forward_2nd_derivatives
                                [point][i][j][l][m] =
                                  tmp[i][j][l][0] * data.covariant[point][m][0];
                              for (unsigned int mr = 1; mr < dim; ++mr)
                                jacobian_pushed_forward_2nd_derivatives[point]
                                                                       [i][j][l]
                                                                       [m] +=
                                  tmp[i][j][l][mr] *
                                  data.covariant[point][m][mr];
                            }
                  }
              }
          }
      }

      /**
       * Update the fourth derivatives of the transformation from unit to real
       * cell, the Jacobian hessian gradients.
       *
       * Skip the computation if possible as indicated by the first argument.
       */
      template <int dim, int spacedim>
      void
      maybe_update_jacobian_3rd_derivatives(
        const CellSimilarity::Similarity                  cell_similarity,
        const typename QProjector<dim>::DataSetDescriptor data_set,
        const typename dealii::MappingFE<dim, spacedim>::InternalData &data,
        std::vector<DerivativeForm<4, dim, spacedim>> &jacobian_3rd_derivatives,
        const unsigned int                             n_q_points)
      {
        const UpdateFlags update_flags = data.update_each;
        if (update_flags & update_jacobian_3rd_derivatives)
          {
            AssertIndexRange(n_q_points, jacobian_3rd_derivatives.size() + 1);

            if (cell_similarity != CellSimilarity::translation)
              {
                for (unsigned int point = 0; point < n_q_points; ++point)
                  {
                    const Tensor<4, dim> *fourth =
                      &data.fourth_derivative(point + data_set, 0);
                    double result[spacedim][dim][dim][dim][dim];
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          for (unsigned int m = 0; m < dim; ++m)
                            for (unsigned int n = 0; n < dim; ++n)
                              result[i][j][l][m][n] =
                                (fourth[0][j][l][m][n] *
                                 data.mapping_support_points[0][i]);
                    for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                      for (unsigned int i = 0; i < spacedim; ++i)
                        for (unsigned int j = 0; j < dim; ++j)
                          for (unsigned int l = 0; l < dim; ++l)
                            for (unsigned int m = 0; m < dim; ++m)
                              for (unsigned int n = 0; n < dim; ++n)
                                result[i][j][l][m][n] +=
                                  (fourth[k][j][l][m][n] *
                                   data.mapping_support_points[k][i]);

                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          for (unsigned int m = 0; m < dim; ++m)
                            for (unsigned int n = 0; n < dim; ++n)
                              jacobian_3rd_derivatives[point][i][j][l][m][n] =
                                result[i][j][l][m][n];
                  }
              }
          }
      }

      /**
       * Update the Hessian gradient of the transformation from unit to real
       * cell, the Jacobian Hessians, pushed forward to the real cell
       * coordinates.
       *
       * Skip the computation if possible as indicated by the first argument.
       */
      template <int dim, int spacedim>
      void
      maybe_update_jacobian_pushed_forward_3rd_derivatives(
        const CellSimilarity::Similarity                  cell_similarity,
        const typename QProjector<dim>::DataSetDescriptor data_set,
        const typename dealii::MappingFE<dim, spacedim>::InternalData &data,
        std::vector<Tensor<5, spacedim>>
                          &jacobian_pushed_forward_3rd_derivatives,
        const unsigned int n_q_points)
      {
        const UpdateFlags update_flags = data.update_each;
        if (update_flags & update_jacobian_pushed_forward_3rd_derivatives)
          {
            AssertIndexRange(n_q_points,
                             jacobian_pushed_forward_3rd_derivatives.size() +
                               1);

            if (cell_similarity != CellSimilarity::translation)
              {
                double tmp[spacedim][spacedim][spacedim][spacedim][spacedim];
                for (unsigned int point = 0; point < n_q_points; ++point)
                  {
                    const Tensor<4, dim> *fourth =
                      &data.fourth_derivative(point + data_set, 0);
                    double result[spacedim][dim][dim][dim][dim];
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          for (unsigned int m = 0; m < dim; ++m)
                            for (unsigned int n = 0; n < dim; ++n)
                              result[i][j][l][m][n] =
                                (fourth[0][j][l][m][n] *
                                 data.mapping_support_points[0][i]);
                    for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                      for (unsigned int i = 0; i < spacedim; ++i)
                        for (unsigned int j = 0; j < dim; ++j)
                          for (unsigned int l = 0; l < dim; ++l)
                            for (unsigned int m = 0; m < dim; ++m)
                              for (unsigned int n = 0; n < dim; ++n)
                                result[i][j][l][m][n] +=
                                  (fourth[k][j][l][m][n] *
                                   data.mapping_support_points[k][i]);

                    // push-forward the j-coordinate
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < spacedim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          for (unsigned int m = 0; m < dim; ++m)
                            for (unsigned int n = 0; n < dim; ++n)
                              {
                                tmp[i][j][l][m][n] =
                                  result[i][0][l][m][n] *
                                  data.covariant[point][j][0];
                                for (unsigned int jr = 1; jr < dim; ++jr)
                                  tmp[i][j][l][m][n] +=
                                    result[i][jr][l][m][n] *
                                    data.covariant[point][j][jr];
                              }

                    // push-forward the l-coordinate
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < spacedim; ++j)
                        for (unsigned int l = 0; l < spacedim; ++l)
                          for (unsigned int m = 0; m < dim; ++m)
                            for (unsigned int n = 0; n < dim; ++n)
                              {
                                jacobian_pushed_forward_3rd_derivatives
                                  [point][i][j][l][m][n] =
                                    tmp[i][j][0][m][n] *
                                    data.covariant[point][l][0];
                                for (unsigned int lr = 1; lr < dim; ++lr)
                                  jacobian_pushed_forward_3rd_derivatives
                                    [point][i][j][l][m][n] +=
                                    tmp[i][j][lr][m][n] *
                                    data.covariant[point][l][lr];
                              }

                    // push-forward the m-coordinate
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < spacedim; ++j)
                        for (unsigned int l = 0; l < spacedim; ++l)
                          for (unsigned int m = 0; m < spacedim; ++m)
                            for (unsigned int n = 0; n < dim; ++n)
                              {
                                tmp[i][j][l][m][n] =
                                  jacobian_pushed_forward_3rd_derivatives
                                    [point][i][j][l][0][n] *
                                  data.covariant[point][m][0];
                                for (unsigned int mr = 1; mr < dim; ++mr)
                                  tmp[i][j][l][m][n] +=
                                    jacobian_pushed_forward_3rd_derivatives
                                      [point][i][j][l][mr][n] *
                                    data.covariant[point][m][mr];
                              }

                    // push-forward the n-coordinate
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < spacedim; ++j)
                        for (unsigned int l = 0; l < spacedim; ++l)
                          for (unsigned int m = 0; m < spacedim; ++m)
                            for (unsigned int n = 0; n < spacedim; ++n)
                              {
                                jacobian_pushed_forward_3rd_derivatives
                                  [point][i][j][l][m][n] =
                                    tmp[i][j][l][m][0] *
                                    data.covariant[point][n][0];
                                for (unsigned int nr = 1; nr < dim; ++nr)
                                  jacobian_pushed_forward_3rd_derivatives
                                    [point][i][j][l][m][n] +=
                                    tmp[i][j][l][m][nr] *
                                    data.covariant[point][n][nr];
                              }
                  }
              }
          }
      }
    } // namespace
  }   // namespace MappingFEImplementation
} // namespace internal



template <int dim, int spacedim>
MappingFE<dim, spacedim>::MappingFE(const FiniteElement<dim, spacedim> &fe)
  : fe(fe.clone())
  , polynomial_degree(fe.tensor_degree())
{
  Assert(polynomial_degree >= 1,
         ExcMessage("It only makes sense to create polynomial mappings "
                    "with a polynomial degree greater or equal to one."));
  Assert(fe.n_components() == 1, ExcNotImplemented());

  Assert(fe.has_support_points(), ExcNotImplemented());

  const auto &mapping_support_points = fe.get_unit_support_points();

  const auto reference_cell = fe.reference_cell();

  const unsigned int n_points          = mapping_support_points.size();
  const unsigned int n_shape_functions = reference_cell.n_vertices();

  this->mapping_support_point_weights =
    Table<2, double>(n_points, n_shape_functions);

  for (unsigned int point = 0; point < n_points; ++point)
    for (unsigned int i = 0; i < n_shape_functions; ++i)
      mapping_support_point_weights(point, i) =
        reference_cell.d_linear_shape_function(mapping_support_points[point],
                                               i);
}



template <int dim, int spacedim>
MappingFE<dim, spacedim>::MappingFE(const MappingFE<dim, spacedim> &mapping)
  : fe(mapping.fe->clone())
  , polynomial_degree(mapping.polynomial_degree)
  , mapping_support_point_weights(mapping.mapping_support_point_weights)
{}



template <int dim, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
MappingFE<dim, spacedim>::clone() const
{
  return std::make_unique<MappingFE<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
unsigned int
MappingFE<dim, spacedim>::get_degree() const
{
  return polynomial_degree;
}



template <int dim, int spacedim>
Point<spacedim>
MappingFE<dim, spacedim>::transform_unit_to_real_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<dim>                                           &p) const
{
  const std::vector<Point<spacedim>> support_points =
    this->compute_mapping_support_points(cell);

  Point<spacedim> mapped_point;

  for (unsigned int i = 0; i < this->fe->n_dofs_per_cell(); ++i)
    mapped_point += support_points[i] * this->fe->shape_value(i, p);

  return mapped_point;
}



template <int dim, int spacedim>
Point<dim>
MappingFE<dim, spacedim>::transform_real_to_unit_cell(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<spacedim>                                      &p) const
{
  const std::vector<Point<spacedim>> support_points =
    this->compute_mapping_support_points(cell);

  const double       eps        = 1.e-12 * cell->diameter();
  const unsigned int loop_limit = 10;

  Point<dim> p_unit;

  unsigned int loop = 0;

  // This loop solves the following problem:
  // grad_F^T residual + (grad_F^T grad_F + grad_F^T hess_F^T dp) dp = 0
  // where the term
  // (grad_F^T hess_F dp) is approximated by (-hess_F * residual)
  // This is basically a second order approximation of Newton method, where the
  // Jacobian is corrected with a higher order term coming from the hessian.
  do
    {
      Point<spacedim> mapped_point;

      // Transpose of the gradient map
      DerivativeForm<1, spacedim, dim> grad_FT;
      DerivativeForm<2, spacedim, dim> hess_FT;

      for (unsigned int i = 0; i < this->fe->n_dofs_per_cell(); ++i)
        {
          mapped_point += support_points[i] * this->fe->shape_value(i, p_unit);
          const auto grad_F_i    = this->fe->shape_grad(i, p_unit);
          const auto hessian_F_i = this->fe->shape_grad_grad(i, p_unit);
          for (unsigned int j = 0; j < dim; ++j)
            {
              grad_FT[j] += grad_F_i[j] * support_points[i];
              for (unsigned int l = 0; l < dim; ++l)
                hess_FT[j][l] += hessian_F_i[j][l] * support_points[i];
            }
        }

      // Residual
      const auto residual = p - mapped_point;
      // Project the residual on the reference coordinate system
      // to compute the error, and to filter components orthogonal to the
      // manifold, and compute a 2nd order correction of the metric tensor
      const auto grad_FT_residual = apply_transformation(grad_FT, residual);

      // Do not invert nor compute the metric if not necessary.
      if (grad_FT_residual.norm() <= eps)
        break;

      // Now compute the (corrected) metric tensor
      Tensor<2, dim> corrected_metric_tensor;
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int l = 0; l < dim; ++l)
          corrected_metric_tensor[j][l] =
            -grad_FT[j] * grad_FT[l] + hess_FT[j][l] * residual;

      // And compute the update
      const auto g_inverse = invert(corrected_metric_tensor);
      p_unit -= Point<dim>(g_inverse * grad_FT_residual);

      ++loop;
    }
  while (loop < loop_limit);

  // Here we check that in the last execution of while the first
  // condition was already wrong, meaning the residual was below
  // eps. Only if the first condition failed, loop will have been
  // increased and tested, and thus have reached the limit.
  AssertThrow(loop < loop_limit,
              (typename Mapping<dim, spacedim>::ExcTransformationFailed()));

  return p_unit;
}



template <int dim, int spacedim>
UpdateFlags
MappingFE<dim, spacedim>::requires_update_flags(const UpdateFlags in) const
{
  // add flags if the respective quantities are necessary to compute
  // what we need. note that some flags appear in both the conditions
  // and in subsequent set operations. this leads to some circular
  // logic. the only way to treat this is to iterate. since there are
  // 5 if-clauses in the loop, it will take at most 5 iterations to
  // converge. do them:
  UpdateFlags out = in;
  for (unsigned int i = 0; i < 5; ++i)
    {
      // The following is a little incorrect:
      // If not applied on a face,
      // update_boundary_forms does not
      // make sense. On the other hand,
      // it is necessary on a
      // face. Currently,
      // update_boundary_forms is simply
      // ignored for the interior of a
      // cell.
      if (out & (update_JxW_values | update_normal_vectors))
        out |= update_boundary_forms;

      if (out & (update_covariant_transformation | update_JxW_values |
                 update_jacobians | update_jacobian_grads |
                 update_boundary_forms | update_normal_vectors))
        out |= update_contravariant_transformation;

      if (out &
          (update_inverse_jacobians | update_jacobian_pushed_forward_grads |
           update_jacobian_pushed_forward_2nd_derivatives |
           update_jacobian_pushed_forward_3rd_derivatives))
        out |= update_covariant_transformation;

      // The contravariant transformation is used in the Piola
      // transformation, which requires the determinant of the Jacobi
      // matrix of the transformation.  Because we have no way of
      // knowing here whether the finite element wants to use the
      // contravariant or the Piola transforms, we add the JxW values
      // to the list of flags to be updated for each cell.
      if (out & update_contravariant_transformation)
        out |= update_volume_elements;

      // the same is true when computing normal vectors: they require
      // the determinant of the Jacobian
      if (out & update_normal_vectors)
        out |= update_volume_elements;
    }

  return out;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingFE<dim, spacedim>::get_data(const UpdateFlags      update_flags,
                                   const Quadrature<dim> &q) const
{
  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase> data_ptr =
    std::make_unique<InternalData>(*this->fe);
  data_ptr->reinit(this->requires_update_flags(update_flags), q);
  return data_ptr;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingFE<dim, spacedim>::get_face_data(
  const UpdateFlags               update_flags,
  const hp::QCollection<dim - 1> &quadrature) const
{
  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase> data_ptr =
    std::make_unique<InternalData>(*this->fe);
  auto &data = dynamic_cast<InternalData &>(*data_ptr);
  data.initialize_face(this->requires_update_flags(update_flags),
                       QProjector<dim>::project_to_all_faces(
                         this->fe->reference_cell(), quadrature),
                       quadrature.max_n_quadrature_points());

  return data_ptr;
}



template <int dim, int spacedim>
std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
MappingFE<dim, spacedim>::get_subface_data(
  const UpdateFlags          update_flags,
  const Quadrature<dim - 1> &quadrature) const
{
  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase> data_ptr =
    std::make_unique<InternalData>(*this->fe);
  auto &data = dynamic_cast<InternalData &>(*data_ptr);
  data.initialize_face(this->requires_update_flags(update_flags),
                       QProjector<dim>::project_to_all_subfaces(
                         this->fe->reference_cell(), quadrature),
                       quadrature.size());

  return data_ptr;
}



template <int dim, int spacedim>
CellSimilarity::Similarity
MappingFE<dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const CellSimilarity::Similarity                            cell_similarity,
  const Quadrature<dim>                                      &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // ensure that the following static_cast is really correct:
  Assert(dynamic_cast<const InternalData *>(&internal_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);

  const unsigned int n_q_points = quadrature.size();

  // recompute the support points of the transformation of this
  // cell. we tried to be clever here in an earlier version of the
  // library by checking whether the cell is the same as the one we
  // had visited last, but it turns out to be difficult to determine
  // that because a cell for the purposes of a mapping is
  // characterized not just by its (triangulation, level, index)
  // triple, but also by the locations of its vertices, the manifold
  // object attached to the cell and all of its bounding faces/edges,
  // etc. to reliably test that the "cell" we are on is, therefore,
  // not easily done
  data.mapping_support_points = this->compute_mapping_support_points(cell);
  data.cell_of_current_support_points = cell;

  // if the order of the mapping is greater than 1, then do not reuse any cell
  // similarity information. This is necessary because the cell similarity
  // value is computed with just cell vertices and does not take into account
  // cell curvature.
  const CellSimilarity::Similarity computed_cell_similarity =
    (polynomial_degree == 1 ? cell_similarity : CellSimilarity::none);

  internal::MappingFEImplementation::maybe_compute_q_points<dim, spacedim>(
    QProjector<dim>::DataSetDescriptor::cell(),
    data,
    output_data.quadrature_points,
    n_q_points);

  internal::MappingFEImplementation::maybe_update_Jacobians<dim, spacedim>(
    computed_cell_similarity,
    QProjector<dim>::DataSetDescriptor::cell(),
    data,
    n_q_points);

  internal::MappingFEImplementation::maybe_update_jacobian_grads<dim, spacedim>(
    computed_cell_similarity,
    QProjector<dim>::DataSetDescriptor::cell(),
    data,
    output_data.jacobian_grads,
    n_q_points);

  internal::MappingFEImplementation::maybe_update_jacobian_pushed_forward_grads<
    dim,
    spacedim>(computed_cell_similarity,
              QProjector<dim>::DataSetDescriptor::cell(),
              data,
              output_data.jacobian_pushed_forward_grads,
              n_q_points);

  internal::MappingFEImplementation::maybe_update_jacobian_2nd_derivatives<
    dim,
    spacedim>(computed_cell_similarity,
              QProjector<dim>::DataSetDescriptor::cell(),
              data,
              output_data.jacobian_2nd_derivatives,
              n_q_points);

  internal::MappingFEImplementation::
    maybe_update_jacobian_pushed_forward_2nd_derivatives<dim, spacedim>(
      computed_cell_similarity,
      QProjector<dim>::DataSetDescriptor::cell(),
      data,
      output_data.jacobian_pushed_forward_2nd_derivatives,
      n_q_points);

  internal::MappingFEImplementation::maybe_update_jacobian_3rd_derivatives<
    dim,
    spacedim>(computed_cell_similarity,
              QProjector<dim>::DataSetDescriptor::cell(),
              data,
              output_data.jacobian_3rd_derivatives,
              n_q_points);

  internal::MappingFEImplementation::
    maybe_update_jacobian_pushed_forward_3rd_derivatives<dim, spacedim>(
      computed_cell_similarity,
      QProjector<dim>::DataSetDescriptor::cell(),
      data,
      output_data.jacobian_pushed_forward_3rd_derivatives,
      n_q_points);

  const UpdateFlags          update_flags = data.update_each;
  const std::vector<double> &weights      = quadrature.get_weights();

  // Multiply quadrature weights by absolute value of Jacobian determinants or
  // the area element g=sqrt(DX^t DX) in case of codim > 0

  if (update_flags & (update_normal_vectors | update_JxW_values))
    {
      AssertDimension(output_data.JxW_values.size(), n_q_points);

      Assert(!(update_flags & update_normal_vectors) ||
               (output_data.normal_vectors.size() == n_q_points),
             ExcDimensionMismatch(output_data.normal_vectors.size(),
                                  n_q_points));


      if (computed_cell_similarity != CellSimilarity::translation)
        for (unsigned int point = 0; point < n_q_points; ++point)
          {
            if (dim == spacedim)
              {
                const double det = data.contravariant[point].determinant();

                // check for distorted cells.

                // TODO: this allows for anisotropies of up to 1e6 in 3d and
                // 1e12 in 2d. might want to find a finer
                // (dimension-independent) criterion
                Assert(det >
                         1e-12 * Utilities::fixed_power<dim>(
                                   cell->diameter() / std::sqrt(double(dim))),
                       (typename Mapping<dim, spacedim>::ExcDistortedMappedCell(
                         cell->center(), det, point)));

                output_data.JxW_values[point] = weights[point] * det;
              }
            // if dim==spacedim, then there is no cell normal to
            // compute. since this is for FEValues (and not FEFaceValues),
            // there are also no face normals to compute
            else // codim>0 case
              {
                Tensor<1, spacedim> DX_t[dim];
                for (unsigned int i = 0; i < spacedim; ++i)
                  for (unsigned int j = 0; j < dim; ++j)
                    DX_t[j][i] = data.contravariant[point][i][j];

                Tensor<2, dim> G; // First fundamental form
                for (unsigned int i = 0; i < dim; ++i)
                  for (unsigned int j = 0; j < dim; ++j)
                    G[i][j] = DX_t[i] * DX_t[j];

                output_data.JxW_values[point] =
                  std::sqrt(determinant(G)) * weights[point];

                if (computed_cell_similarity ==
                    CellSimilarity::inverted_translation)
                  {
                    // we only need to flip the normal
                    if (update_flags & update_normal_vectors)
                      output_data.normal_vectors[point] *= -1.;
                  }
                else
                  {
                    if (update_flags & update_normal_vectors)
                      {
                        Assert(spacedim == dim + 1,
                               ExcMessage(
                                 "There is no (unique) cell normal for " +
                                 Utilities::int_to_string(dim) +
                                 "-dimensional cells in " +
                                 Utilities::int_to_string(spacedim) +
                                 "-dimensional space. This only works if the "
                                 "space dimension is one greater than the "
                                 "dimensionality of the mesh cells."));

                        if (dim == 1)
                          output_data.normal_vectors[point] =
                            cross_product_2d(-DX_t[0]);
                        else // dim == 2
                          output_data.normal_vectors[point] =
                            cross_product_3d(DX_t[0], DX_t[1]);

                        output_data.normal_vectors[point] /=
                          output_data.normal_vectors[point].norm();

                        if (cell->direction_flag() == false)
                          output_data.normal_vectors[point] *= -1.;
                      }
                  }
              } // codim>0 case
          }
    }



  // copy values from InternalData to vector given by reference
  if (update_flags & update_jacobians)
    {
      AssertDimension(output_data.jacobians.size(), n_q_points);
      if (computed_cell_similarity != CellSimilarity::translation)
        for (unsigned int point = 0; point < n_q_points; ++point)
          output_data.jacobians[point] = data.contravariant[point];
    }

  // copy values from InternalData to vector given by reference
  if (update_flags & update_inverse_jacobians)
    {
      AssertDimension(output_data.inverse_jacobians.size(), n_q_points);
      if (computed_cell_similarity != CellSimilarity::translation)
        for (unsigned int point = 0; point < n_q_points; ++point)
          output_data.inverse_jacobians[point] =
            data.covariant[point].transpose();
    }

  return computed_cell_similarity;
}



namespace internal
{
  namespace MappingFEImplementation
  {
    namespace
    {
      /**
       * Depending on what information is called for in the update flags of
       * the
       * @p data object, compute the various pieces of information that is
       * required by the fill_fe_face_values() and fill_fe_subface_values()
       * functions. This function simply unifies the work that would be done
       * by those two functions.
       *
       * The resulting data is put into the @p output_data argument.
       */
      template <int dim, int spacedim>
      void
      maybe_compute_face_data(
        const dealii::MappingFE<dim, spacedim> &mapping,
        const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                                                         &cell,
        const unsigned int                                face_no,
        const unsigned int                                subface_no,
        const unsigned int                                n_q_points,
        const typename QProjector<dim>::DataSetDescriptor data_set,
        const typename dealii::MappingFE<dim, spacedim>::InternalData &data,
        internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
          &output_data)
      {
        const UpdateFlags update_flags = data.update_each;

        if (update_flags &
            (update_boundary_forms | update_normal_vectors | update_jacobians |
             update_JxW_values | update_inverse_jacobians))
          {
            if (update_flags & update_boundary_forms)
              AssertIndexRange(n_q_points,
                               output_data.boundary_forms.size() + 1);
            if (update_flags & update_normal_vectors)
              AssertIndexRange(n_q_points,
                               output_data.normal_vectors.size() + 1);
            if (update_flags & update_JxW_values)
              AssertIndexRange(n_q_points, output_data.JxW_values.size() + 1);

            Assert(data.aux.size() + 1 >= dim, ExcInternalError());

            // first compute some common data that is used for evaluating
            // all of the flags below

            // map the unit tangentials to the real cell. checking for
            // d!=dim-1 eliminates compiler warnings regarding unsigned int
            // expressions < 0.
            for (unsigned int d = 0; d != dim - 1; ++d)
              {
                Assert(face_no + cell->n_faces() * d <
                         data.unit_tangentials.size(),
                       ExcInternalError());
                Assert(
                  data.aux[d].size() <=
                    data.unit_tangentials[face_no + cell->n_faces() * d].size(),
                  ExcInternalError());

                mapping.transform(
                  make_array_view(
                    data.unit_tangentials[face_no + cell->n_faces() * d]),
                  mapping_contravariant,
                  data,
                  make_array_view(data.aux[d]));
              }

            if (update_flags & update_boundary_forms)
              {
                // if dim==spacedim, we can use the unit tangentials to
                // compute the boundary form by simply taking the cross
                // product
                if (dim == spacedim)
                  {
                    for (unsigned int i = 0; i < n_q_points; ++i)
                      switch (dim)
                        {
                          case 1:
                            // in 1d, we don't have access to any of the
                            // data.aux fields (because it has only dim-1
                            // components), but we can still compute the
                            // boundary form by simply looking at the number
                            // of the face
                            output_data.boundary_forms[i][0] =
                              (face_no == 0 ? -1 : +1);
                            break;
                          case 2:
                            output_data.boundary_forms[i] =
                              cross_product_2d(data.aux[0][i]);
                            break;
                          case 3:
                            output_data.boundary_forms[i] =
                              cross_product_3d(data.aux[0][i], data.aux[1][i]);
                            break;
                          default:
                            DEAL_II_NOT_IMPLEMENTED();
                        }
                  }
                else //(dim < spacedim)
                  {
                    // in the codim-one case, the boundary form results from
                    // the cross product of all the face tangential vectors
                    // and the cell normal vector
                    //
                    // to compute the cell normal, use the same method used in
                    // fill_fe_values for cells above
                    AssertIndexRange(n_q_points, data.contravariant.size() + 1);

                    for (unsigned int point = 0; point < n_q_points; ++point)
                      {
                        if (dim == 1)
                          {
                            // J is a tangent vector
                            output_data.boundary_forms[point] =
                              data.contravariant[point].transpose()[0];
                            output_data.boundary_forms[point] /=
                              (face_no == 0 ? -1. : +1.) *
                              output_data.boundary_forms[point].norm();
                          }

                        if (dim == 2)
                          {
                            const DerivativeForm<1, spacedim, dim> DX_t =
                              data.contravariant[point].transpose();

                            Tensor<1, spacedim> cell_normal =
                              cross_product_3d(DX_t[0], DX_t[1]);
                            cell_normal /= cell_normal.norm();

                            // then compute the face normal from the face
                            // tangent and the cell normal:
                            output_data.boundary_forms[point] =
                              cross_product_3d(data.aux[0][point], cell_normal);
                          }
                      }
                  }
              }

            if (update_flags & update_JxW_values)
              for (unsigned int i = 0; i < n_q_points; ++i)
                {
                  output_data.JxW_values[i] =
                    output_data.boundary_forms[i].norm() *
                    data.quadrature_weights[i + data_set];

                  if (subface_no != numbers::invalid_unsigned_int)
                    {
                      if (dim == 2)
                        {
                          const double area_ratio =
                            1. / cell->reference_cell()
                                   .face_reference_cell(face_no)
                                   .n_isotropic_children();
                          output_data.JxW_values[i] *= area_ratio;
                        }
                      else
                        DEAL_II_NOT_IMPLEMENTED();
                    }
                }

            if (update_flags & update_normal_vectors)
              for (unsigned int i = 0; i < n_q_points; ++i)
                output_data.normal_vectors[i] =
                  Point<spacedim>(output_data.boundary_forms[i] /
                                  output_data.boundary_forms[i].norm());

            if (update_flags & update_jacobians)
              for (unsigned int point = 0; point < n_q_points; ++point)
                output_data.jacobians[point] = data.contravariant[point];

            if (update_flags & update_inverse_jacobians)
              for (unsigned int point = 0; point < n_q_points; ++point)
                output_data.inverse_jacobians[point] =
                  data.covariant[point].transpose();
          }
      }


      /**
       * Do the work of MappingFE::fill_fe_face_values() and
       * MappingFE::fill_fe_subface_values() in a generic way,
       * using the 'data_set' to differentiate whether we will
       * work on a face (and if so, which one) or subface.
       */
      template <int dim, int spacedim>
      void
      do_fill_fe_face_values(
        const dealii::MappingFE<dim, spacedim> &mapping,
        const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                                                         &cell,
        const unsigned int                                face_no,
        const unsigned int                                subface_no,
        const typename QProjector<dim>::DataSetDescriptor data_set,
        const Quadrature<dim - 1>                        &quadrature,
        const typename dealii::MappingFE<dim, spacedim>::InternalData &data,
        internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
          &output_data)
      {
        const unsigned int n_q_points = quadrature.size();

        maybe_compute_q_points<dim, spacedim>(data_set,
                                              data,
                                              output_data.quadrature_points,
                                              n_q_points);
        maybe_update_Jacobians<dim, spacedim>(CellSimilarity::none,
                                              data_set,
                                              data,
                                              n_q_points);
        maybe_update_jacobian_grads<dim, spacedim>(CellSimilarity::none,
                                                   data_set,
                                                   data,
                                                   output_data.jacobian_grads,
                                                   n_q_points);
        maybe_update_jacobian_pushed_forward_grads<dim, spacedim>(
          CellSimilarity::none,
          data_set,
          data,
          output_data.jacobian_pushed_forward_grads,
          n_q_points);
        maybe_update_jacobian_2nd_derivatives<dim, spacedim>(
          CellSimilarity::none,
          data_set,
          data,
          output_data.jacobian_2nd_derivatives,
          n_q_points);
        maybe_update_jacobian_pushed_forward_2nd_derivatives<dim, spacedim>(
          CellSimilarity::none,
          data_set,
          data,
          output_data.jacobian_pushed_forward_2nd_derivatives,
          n_q_points);
        maybe_update_jacobian_3rd_derivatives<dim, spacedim>(
          CellSimilarity::none,
          data_set,
          data,
          output_data.jacobian_3rd_derivatives,
          n_q_points);
        maybe_update_jacobian_pushed_forward_3rd_derivatives<dim, spacedim>(
          CellSimilarity::none,
          data_set,
          data,
          output_data.jacobian_pushed_forward_3rd_derivatives,
          n_q_points);

        maybe_compute_face_data(mapping,
                                cell,
                                face_no,
                                subface_no,
                                n_q_points,
                                data_set,
                                data,
                                output_data);
      }
    } // namespace
  }   // namespace MappingFEImplementation
} // namespace internal



template <int dim, int spacedim>
void
MappingFE<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const hp::QCollection<dim - 1>                             &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // ensure that the following cast is really correct:
  Assert((dynamic_cast<const InternalData *>(&internal_data) != nullptr),
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);

  // if necessary, recompute the support points of the transformation of this
  // cell (note that we need to first check the triangulation pointer, since
  // otherwise the second test might trigger an exception if the
  // triangulations are not the same)
  if ((data.mapping_support_points.empty()) ||
      (&cell->get_triangulation() !=
       &data.cell_of_current_support_points->get_triangulation()) ||
      (cell != data.cell_of_current_support_points))
    {
      data.mapping_support_points = this->compute_mapping_support_points(cell);
      data.cell_of_current_support_points = cell;
    }

  internal::MappingFEImplementation::do_fill_fe_face_values(
    *this,
    cell,
    face_no,
    numbers::invalid_unsigned_int,
    QProjector<dim>::DataSetDescriptor::face(this->fe->reference_cell(),
                                             face_no,
                                             cell->combined_face_orientation(
                                               face_no),
                                             quadrature),
    quadrature[quadrature.size() == 1 ? 0 : face_no],
    data,
    output_data);
}



template <int dim, int spacedim>
void
MappingFE<dim, spacedim>::fill_fe_subface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const unsigned int                                          subface_no,
  const Quadrature<dim - 1>                                  &quadrature,
  const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &output_data) const
{
  // ensure that the following cast is really correct:
  Assert((dynamic_cast<const InternalData *>(&internal_data) != nullptr),
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);

  // if necessary, recompute the support points of the transformation of this
  // cell (note that we need to first check the triangulation pointer, since
  // otherwise the second test might trigger an exception if the
  // triangulations are not the same)
  if ((data.mapping_support_points.empty()) ||
      (&cell->get_triangulation() !=
       &data.cell_of_current_support_points->get_triangulation()) ||
      (cell != data.cell_of_current_support_points))
    {
      data.mapping_support_points = this->compute_mapping_support_points(cell);
      data.cell_of_current_support_points = cell;
    }

  internal::MappingFEImplementation::do_fill_fe_face_values(
    *this,
    cell,
    face_no,
    subface_no,
    QProjector<dim>::DataSetDescriptor::subface(this->fe->reference_cell(),
                                                face_no,
                                                subface_no,
                                                cell->combined_face_orientation(
                                                  face_no),
                                                quadrature.size(),
                                                cell->subface_case(face_no)),
    quadrature,
    data,
    output_data);
}



namespace internal
{
  namespace MappingFEImplementation
  {
    namespace
    {
      template <int dim, int spacedim, int rank>
      void
      transform_fields(
        const ArrayView<const Tensor<rank, dim>>                &input,
        const MappingKind                                        mapping_kind,
        const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
        const ArrayView<Tensor<rank, spacedim>>                 &output)
      {
        // In the case of wedges and pyramids, faces might have different
        // numbers of quadrature points on each face with the result
        // that input and output have different sizes, since input has
        // the correct size but the size of output is the maximum of
        // all possible sizes.
        AssertIndexRange(input.size(), output.size() + 1);

        Assert(
          (dynamic_cast<
             const typename dealii::MappingFE<dim, spacedim>::InternalData *>(
             &mapping_data) != nullptr),
          ExcInternalError());
        const typename dealii::MappingFE<dim, spacedim>::InternalData &data =
          static_cast<
            const typename dealii::MappingFE<dim, spacedim>::InternalData &>(
            mapping_data);

        switch (mapping_kind)
          {
            case mapping_contravariant:
              {
                Assert(
                  data.update_each & update_contravariant_transformation,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_contravariant_transformation"));

                for (unsigned int i = 0; i < input.size(); ++i)
                  output[i] =
                    apply_transformation(data.contravariant[i], input[i]);

                return;
              }

            case mapping_piola:
              {
                Assert(
                  data.update_each & update_contravariant_transformation,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_contravariant_transformation"));
                Assert(
                  data.update_each & update_volume_elements,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_volume_elements"));
                Assert(rank == 1, ExcMessage("Only for rank 1"));
                if (rank != 1)
                  return;

                for (unsigned int i = 0; i < input.size(); ++i)
                  {
                    output[i] =
                      apply_transformation(data.contravariant[i], input[i]);
                    output[i] /= data.volume_elements[i];
                  }
                return;
              }
            // We still allow this operation as in the
            // reference cell Derivatives are Tensor
            // rather than DerivativeForm
            case mapping_covariant:
              {
                Assert(
                  data.update_each & update_contravariant_transformation,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_covariant_transformation"));

                for (unsigned int i = 0; i < input.size(); ++i)
                  output[i] = apply_transformation(data.covariant[i], input[i]);

                return;
              }

            default:
              DEAL_II_NOT_IMPLEMENTED();
          }
      }


      template <int dim, int spacedim, int rank>
      void
      transform_gradients(
        const ArrayView<const Tensor<rank, dim>>                &input,
        const MappingKind                                        mapping_kind,
        const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
        const ArrayView<Tensor<rank, spacedim>>                 &output)
      {
        AssertDimension(input.size(), output.size());
        Assert(
          (dynamic_cast<
             const typename dealii::MappingFE<dim, spacedim>::InternalData *>(
             &mapping_data) != nullptr),
          ExcInternalError());
        const typename dealii::MappingFE<dim, spacedim>::InternalData &data =
          static_cast<
            const typename dealii::MappingFE<dim, spacedim>::InternalData &>(
            mapping_data);

        switch (mapping_kind)
          {
            case mapping_contravariant_gradient:
              {
                Assert(
                  data.update_each & update_covariant_transformation,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_covariant_transformation"));
                Assert(
                  data.update_each & update_contravariant_transformation,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_contravariant_transformation"));
                Assert(rank == 2, ExcMessage("Only for rank 2"));

                for (unsigned int i = 0; i < output.size(); ++i)
                  {
                    const DerivativeForm<1, spacedim, dim> A =
                      apply_transformation(data.contravariant[i],
                                           transpose(input[i]));
                    output[i] =
                      apply_transformation(data.covariant[i], A.transpose());
                  }

                return;
              }

            case mapping_covariant_gradient:
              {
                Assert(
                  data.update_each & update_covariant_transformation,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_covariant_transformation"));
                Assert(rank == 2, ExcMessage("Only for rank 2"));

                for (unsigned int i = 0; i < output.size(); ++i)
                  {
                    const DerivativeForm<1, spacedim, dim> A =
                      apply_transformation(data.covariant[i],
                                           transpose(input[i]));
                    output[i] =
                      apply_transformation(data.covariant[i], A.transpose());
                  }

                return;
              }

            case mapping_piola_gradient:
              {
                Assert(
                  data.update_each & update_covariant_transformation,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_covariant_transformation"));
                Assert(
                  data.update_each & update_contravariant_transformation,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_contravariant_transformation"));
                Assert(
                  data.update_each & update_volume_elements,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_volume_elements"));
                Assert(rank == 2, ExcMessage("Only for rank 2"));

                for (unsigned int i = 0; i < output.size(); ++i)
                  {
                    const DerivativeForm<1, spacedim, dim> A =
                      apply_transformation(data.covariant[i], input[i]);
                    const Tensor<2, spacedim> T =
                      apply_transformation(data.contravariant[i],
                                           A.transpose());

                    output[i] = transpose(T);
                    output[i] /= data.volume_elements[i];
                  }

                return;
              }

            default:
              DEAL_II_NOT_IMPLEMENTED();
          }
      }



      template <int dim, int spacedim>
      void
      transform_hessians(
        const ArrayView<const Tensor<3, dim>>                   &input,
        const MappingKind                                        mapping_kind,
        const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
        const ArrayView<Tensor<3, spacedim>>                    &output)
      {
        AssertDimension(input.size(), output.size());
        Assert(
          (dynamic_cast<
             const typename dealii::MappingFE<dim, spacedim>::InternalData *>(
             &mapping_data) != nullptr),
          ExcInternalError());
        const typename dealii::MappingFE<dim, spacedim>::InternalData &data =
          static_cast<
            const typename dealii::MappingFE<dim, spacedim>::InternalData &>(
            mapping_data);

        switch (mapping_kind)
          {
            case mapping_contravariant_hessian:
              {
                Assert(
                  data.update_each & update_covariant_transformation,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_covariant_transformation"));
                Assert(
                  data.update_each & update_contravariant_transformation,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_contravariant_transformation"));

                for (unsigned int q = 0; q < output.size(); ++q)
                  output[q] =
                    internal::apply_contravariant_hessian(data.covariant[q],
                                                          data.contravariant[q],
                                                          input[q]);

                return;
              }

            case mapping_covariant_hessian:
              {
                Assert(
                  data.update_each & update_covariant_transformation,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_covariant_transformation"));

                for (unsigned int q = 0; q < output.size(); ++q)
                  output[q] =
                    internal::apply_covariant_hessian(data.covariant[q],
                                                      input[q]);

                return;
              }

            case mapping_piola_hessian:
              {
                Assert(
                  data.update_each & update_covariant_transformation,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_covariant_transformation"));
                Assert(
                  data.update_each & update_contravariant_transformation,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_contravariant_transformation"));
                Assert(
                  data.update_each & update_volume_elements,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_volume_elements"));

                for (unsigned int q = 0; q < output.size(); ++q)
                  output[q] =
                    internal::apply_piola_hessian(data.covariant[q],
                                                  data.contravariant[q],
                                                  data.volume_elements[q],
                                                  input[q]);

                return;
              }

            default:
              DEAL_II_NOT_IMPLEMENTED();
          }
      }



      template <int dim, int spacedim, int rank>
      void
      transform_differential_forms(
        const ArrayView<const DerivativeForm<rank, dim, spacedim>> &input,
        const MappingKind                                        mapping_kind,
        const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
        const ArrayView<Tensor<rank + 1, spacedim>>             &output)
      {
        AssertDimension(input.size(), output.size());
        Assert(
          (dynamic_cast<
             const typename dealii::MappingFE<dim, spacedim>::InternalData *>(
             &mapping_data) != nullptr),
          ExcInternalError());
        const typename dealii::MappingFE<dim, spacedim>::InternalData &data =
          static_cast<
            const typename dealii::MappingFE<dim, spacedim>::InternalData &>(
            mapping_data);

        switch (mapping_kind)
          {
            case mapping_covariant:
              {
                Assert(
                  data.update_each & update_contravariant_transformation,
                  typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                    "update_covariant_transformation"));

                for (unsigned int i = 0; i < output.size(); ++i)
                  output[i] = apply_transformation(data.covariant[i], input[i]);

                return;
              }
            default:
              DEAL_II_NOT_IMPLEMENTED();
          }
      }
    } // namespace
  }   // namespace MappingFEImplementation
} // namespace internal



template <int dim, int spacedim>
void
MappingFE<dim, spacedim>::transform(
  const ArrayView<const Tensor<1, dim>>                   &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<1, spacedim>>                    &output) const
{
  internal::MappingFEImplementation::transform_fields(input,
                                                      mapping_kind,
                                                      mapping_data,
                                                      output);
}



template <int dim, int spacedim>
void
MappingFE<dim, spacedim>::transform(
  const ArrayView<const DerivativeForm<1, dim, spacedim>> &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<2, spacedim>>                    &output) const
{
  internal::MappingFEImplementation::transform_differential_forms(input,
                                                                  mapping_kind,
                                                                  mapping_data,
                                                                  output);
}



template <int dim, int spacedim>
void
MappingFE<dim, spacedim>::transform(
  const ArrayView<const Tensor<2, dim>>                   &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<2, spacedim>>                    &output) const
{
  switch (mapping_kind)
    {
      case mapping_contravariant:
        internal::MappingFEImplementation::transform_fields(input,
                                                            mapping_kind,
                                                            mapping_data,
                                                            output);
        return;

      case mapping_piola_gradient:
      case mapping_contravariant_gradient:
      case mapping_covariant_gradient:
        internal::MappingFEImplementation::transform_gradients(input,
                                                               mapping_kind,
                                                               mapping_data,
                                                               output);
        return;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim, int spacedim>
void
MappingFE<dim, spacedim>::transform(
  const ArrayView<const DerivativeForm<2, dim, spacedim>> &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<3, spacedim>>                    &output) const
{
  AssertDimension(input.size(), output.size());
  Assert(dynamic_cast<const InternalData *>(&mapping_data) != nullptr,
         ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(mapping_data);

  switch (mapping_kind)
    {
      case mapping_covariant_gradient:
        {
          Assert(data.update_each & update_covariant_transformation,
                 typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                   "update_covariant_transformation"));

          for (unsigned int q = 0; q < output.size(); ++q)
            output[q] =
              internal::apply_covariant_gradient(data.covariant[q], input[q]);

          return;
        }

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim, int spacedim>
void
MappingFE<dim, spacedim>::transform(
  const ArrayView<const Tensor<3, dim>>                   &input,
  const MappingKind                                        mapping_kind,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  const ArrayView<Tensor<3, spacedim>>                    &output) const
{
  switch (mapping_kind)
    {
      case mapping_piola_hessian:
      case mapping_contravariant_hessian:
      case mapping_covariant_hessian:
        internal::MappingFEImplementation::transform_hessians(input,
                                                              mapping_kind,
                                                              mapping_data,
                                                              output);
        return;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



namespace
{
  template <int spacedim>
  bool
  check_all_manifold_ids_identical(
    const TriaIterator<CellAccessor<1, spacedim>> &)
  {
    return true;
  }



  template <int spacedim>
  bool
  check_all_manifold_ids_identical(
    const TriaIterator<CellAccessor<2, spacedim>> &cell)
  {
    const auto m_id = cell->manifold_id();

    for (const auto f : cell->face_indices())
      if (m_id != cell->face(f)->manifold_id())
        return false;

    return true;
  }



  template <int spacedim>
  bool
  check_all_manifold_ids_identical(
    const TriaIterator<CellAccessor<3, spacedim>> &cell)
  {
    const auto m_id = cell->manifold_id();

    for (const auto f : cell->face_indices())
      if (m_id != cell->face(f)->manifold_id())
        return false;

    for (const auto l : cell->line_indices())
      if (m_id != cell->line(l)->manifold_id())
        return false;

    return true;
  }
} // namespace



template <int dim, int spacedim>
std::vector<Point<spacedim>>
MappingFE<dim, spacedim>::compute_mapping_support_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
{
  std::vector<Point<spacedim>> mapping_support_points(
    fe->get_unit_support_points().size());

  std::vector<Point<spacedim>> vertices(cell->n_vertices());
  for (const unsigned int i : cell->vertex_indices())
    vertices[i] = cell->vertex(i);

  if (check_all_manifold_ids_identical(cell))
    {
      // version 1) all geometric entities have the same manifold:

      cell->get_manifold().get_new_points(vertices,
                                          mapping_support_point_weights,
                                          mapping_support_points);
    }
  else
    {
      // version 1) geometric entities have different manifold

      // helper function to compute mapped points on subentities
      // note[PM]: this function currently only uses the vertices
      // of cells to create new points on subentities; however,
      // one should use all bounding points to create new points
      // as in the case of MappingQ.
      const auto process = [&](const auto    &manifold,
                               const auto    &indices,
                               const unsigned n_points) {
        if ((indices.size() == 0) || (n_points == 0))
          return;

        const unsigned int n_shape_functions =
          this->fe->reference_cell().n_vertices();

        Table<2, double> mapping_support_point_weights_local(n_points,
                                                             n_shape_functions);
        std::vector<Point<spacedim>> mapping_support_points_local(n_points);

        for (unsigned int p = 0; p < n_points; ++p)
          for (unsigned int i = 0; i < n_shape_functions; ++i)
            mapping_support_point_weights_local(p, i) =
              mapping_support_point_weights(
                indices[p + (indices.size() - n_points)], i);

        manifold.get_new_points(vertices,
                                mapping_support_point_weights_local,
                                mapping_support_points_local);

        for (unsigned int p = 0; p < n_points; ++p)
          mapping_support_points[indices[p + (indices.size() - n_points)]] =
            mapping_support_points_local[p];
      };

      // create dummy DoFHandler to extract indices on subobjects
      const auto                  &fe = *this->fe;
      Triangulation<dim, spacedim> tria;
      GridGenerator::reference_cell(tria, fe.reference_cell());
      DoFHandler<dim, spacedim> dof_handler(tria);
      dof_handler.distribute_dofs(fe);
      const auto &cell_ref = dof_handler.begin_active();

      std::vector<types::global_dof_index> indices;

      // add vertices
      for (const unsigned int i : cell->vertex_indices())
        mapping_support_points[i] = cell->vertex(i);

      // process and add line support points
      for (unsigned int l = 0; l < cell_ref->n_lines(); ++l)
        {
          const auto accessor = cell_ref->line(l);
          indices.resize(fe.n_dofs_per_line() + 2 * fe.n_dofs_per_vertex());
          accessor->get_dof_indices(indices);
          process(cell->line(l)->get_manifold(), indices, fe.n_dofs_per_line());
        }

      // process and add face support points
      if constexpr (dim >= 3)
        {
          for (unsigned int f = 0; f < cell_ref->n_faces(); ++f)
            {
              const auto accessor = cell_ref->face(f);
              indices.resize(fe.n_dofs_per_face());
              accessor->get_dof_indices(indices);
              process(cell->face(f)->get_manifold(),
                      indices,
                      fe.n_dofs_per_quad());
            }
        }

      // process and add volume support points
      if constexpr (dim >= 2)
        {
          indices.resize(fe.n_dofs_per_cell());
          cell_ref->get_dof_indices(indices);
          process(cell->get_manifold(),
                  indices,
                  (dim == 2) ? fe.n_dofs_per_quad() : fe.n_dofs_per_hex());
        }
    }

  return mapping_support_points;
}



template <int dim, int spacedim>
BoundingBox<spacedim>
MappingFE<dim, spacedim>::get_bounding_box(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
{
  return BoundingBox<spacedim>(this->compute_mapping_support_points(cell));
}



template <int dim, int spacedim>
bool
MappingFE<dim, spacedim>::is_compatible_with(
  const ReferenceCell &reference_cell) const
{
  Assert(dim == reference_cell.get_dimension(),
         ExcMessage("The dimension of your mapping (" +
                    Utilities::to_string(dim) +
                    ") and the reference cell cell_type (" +
                    Utilities::to_string(reference_cell.get_dimension()) +
                    " ) do not agree."));

  return fe->reference_cell() == reference_cell;
}



//--------------------------- Explicit instantiations -----------------------
#include "fe/mapping_fe.inst"


DEAL_II_NAMESPACE_CLOSE
