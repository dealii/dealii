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


#ifndef dealii_matrix_free_evaluation_kernels_h
#define dealii_matrix_free_evaluation_kernels_h

#include <deal.II/base/config.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/evaluation_kernels_common.h>
#include <deal.II/matrix_free/fe_evaluation_data.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  // Select evaluator type from element shape function type
  template <MatrixFreeFunctions::ElementType element, bool is_long>
  struct EvaluatorSelector
  {};

  template <bool is_long>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_general, is_long>
  {
    static const EvaluatorVariant variant = evaluate_general;
  };

  template <>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric, false>
  {
    static const EvaluatorVariant variant = evaluate_symmetric;
  };

  template <>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric, true>
  {
    static const EvaluatorVariant variant = evaluate_evenodd;
  };

  template <bool is_long>
  struct EvaluatorSelector<MatrixFreeFunctions::truncated_tensor, is_long>
  {
    static const EvaluatorVariant variant = evaluate_general;
  };

  template <>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric_plus_dg0,
                           false>
  {
    static const EvaluatorVariant variant = evaluate_general;
  };

  template <>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric_plus_dg0, true>
  {
    static const EvaluatorVariant variant = evaluate_evenodd;
  };

  template <bool is_long>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric_collocation,
                           is_long>
  {
    static const EvaluatorVariant variant = evaluate_evenodd;
  };



  /**
   * This struct performs the evaluation of function values and gradients for
   * tensor-product finite elements. The operation is used for both the
   * symmetric and non-symmetric case, which use different apply functions
   * 'values', 'gradients' in the individual coordinate directions. The apply
   * functions for values are provided through one of the template classes
   * EvaluatorTensorProduct which in turn are selected from the
   * MatrixFreeFunctions::ElementType template argument.
   *
   * There are two specialized implementation classes
   * FEEvaluationImplCollocation (for Gauss-Lobatto elements where the nodal
   * points and the quadrature points coincide and the 'values' operation is
   * identity) and FEEvaluationImplTransformToCollocation (which can be
   * transformed to a collocation space and can then use the identity in these
   * spaces), which both allow for shorter code.
   *
   * @note Hessians of the solution are handled in the general
   * FEEvaluationImplSelector struct below, because they can be implemented
   * with the only two code paths for all supported cases, including the
   * specialized cases below.
   */
  template <MatrixFreeFunctions::ElementType type,
            int                              dim,
            int                              fe_degree,
            int                              n_q_points_1d,
            typename Number>
  struct FEEvaluationImpl
  {
    static const EvaluatorVariant variant =
      EvaluatorSelector<type, (fe_degree + n_q_points_1d > 4)>::variant;
    using Number2 =
      typename FEEvaluationData<dim, Number, false>::shape_info_number_type;

    using Eval = EvaluatorTensorProduct<variant,
                                        dim,
                                        fe_degree + 1,
                                        n_q_points_1d,
                                        Number,
                                        Number2>;

    static void
    evaluate(const unsigned int                     n_components,
             const EvaluationFlags::EvaluationFlags evaluation_flag,
             const Number                          *values_dofs_actual,
             FEEvaluationData<dim, Number, false>  &fe_eval);

    static void
    integrate(const unsigned int                     n_components,
              const EvaluationFlags::EvaluationFlags integration_flag,
              Number                                *values_dofs_actual,
              FEEvaluationData<dim, Number, false>  &fe_eval,
              const bool                             add_into_values_array);

    static Eval
    create_evaluator_tensor_product(
      const MatrixFreeFunctions::UnivariateShapeData<Number2>
        *univariate_shape_data)
    {
      if (variant == evaluate_evenodd)
        return Eval(univariate_shape_data->shape_values_eo,
                    univariate_shape_data->shape_gradients_eo,
                    univariate_shape_data->shape_hessians_eo,
                    univariate_shape_data->fe_degree + 1,
                    univariate_shape_data->n_q_points_1d);
      else
        return Eval(univariate_shape_data->shape_values,
                    univariate_shape_data->shape_gradients,
                    univariate_shape_data->shape_hessians,
                    univariate_shape_data->fe_degree + 1,
                    univariate_shape_data->n_q_points_1d);
    }
  };



  /**
   * Specialization for MatrixFreeFunctions::tensor_none, which cannot use the
   * sum-factorization kernels.
   */
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct FEEvaluationImpl<MatrixFreeFunctions::tensor_none,
                          dim,
                          fe_degree,
                          n_q_points_1d,
                          Number>
  {
    static void
    evaluate(const unsigned int                     n_components,
             const EvaluationFlags::EvaluationFlags evaluation_flag,
             const Number                          *values_dofs_actual,
             FEEvaluationData<dim, Number, false>  &fe_eval);

    static void
    integrate(const unsigned int                     n_components,
              const EvaluationFlags::EvaluationFlags integration_flag,
              Number                                *values_dofs_actual,
              FEEvaluationData<dim, Number, false>  &fe_eval,
              const bool                             add_into_values_array);
  };



  template <MatrixFreeFunctions::ElementType type,
            int                              dim,
            int                              fe_degree,
            int                              n_q_points_1d,
            typename Number>
  inline void
  FEEvaluationImpl<type, dim, fe_degree, n_q_points_1d, Number>::evaluate(
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags evaluation_flag,
    const Number                          *values_dofs_actual,
    FEEvaluationData<dim, Number, false>  &fe_eval)
  {
    if (evaluation_flag == EvaluationFlags::nothing)
      return;

    std::array<const MatrixFreeFunctions::UnivariateShapeData<Number2> *, 3>
      univariate_shape_data;

    const auto &shape_data = fe_eval.get_shape_info().data;

    univariate_shape_data.fill(&shape_data.front());

    if (shape_data.size() == dim)
      for (int i = 1; i < dim; ++i)
        univariate_shape_data[i] = &shape_data[i];

    Eval eval0 = create_evaluator_tensor_product(univariate_shape_data[0]);
    Eval eval1 = create_evaluator_tensor_product(univariate_shape_data[1]);
    Eval eval2 = create_evaluator_tensor_product(univariate_shape_data[2]);

    const unsigned int temp_size =
      Eval::n_rows_of_product == numbers::invalid_unsigned_int ?
        0 :
        (Eval::n_rows_of_product > Eval::n_columns_of_product ?
           Eval::n_rows_of_product :
           Eval::n_columns_of_product);
    Number *temp1 = fe_eval.get_scratch_data().begin();
    Number *temp2;
    if (temp_size == 0)
      {
        temp2 = temp1 + std::max(Utilities::fixed_power<dim>(
                                   shape_data.front().fe_degree + 1),
                                 Utilities::fixed_power<dim>(
                                   shape_data.front().n_q_points_1d));
      }
    else
      {
        temp2 = temp1 + temp_size;
      }

    const std::size_t n_q_points = temp_size == 0 ?
                                     fe_eval.get_shape_info().n_q_points :
                                     Eval::n_columns_of_product;
    const std::size_t dofs_per_comp =
      (type == MatrixFreeFunctions::truncated_tensor) ?
        Utilities::pow(shape_data.front().fe_degree + 1, dim) :
        fe_eval.get_shape_info().dofs_per_component_on_cell;
    const Number *values_dofs =
      (type == MatrixFreeFunctions::truncated_tensor) ?
        temp1 + 2 * (std::max<std::size_t>(
                      fe_eval.get_shape_info().dofs_per_component_on_cell,
                      n_q_points)) :
        values_dofs_actual;

    if (type == MatrixFreeFunctions::truncated_tensor)
      embed_truncated_into_full_tensor_product<dim, fe_degree>(
        n_components,
        const_cast<Number *>(values_dofs),
        values_dofs_actual,
        fe_eval);

    Number *values_quad    = fe_eval.begin_values();
    Number *gradients_quad = fe_eval.begin_gradients();

    switch (dim)
      {
        case 1:
          for (unsigned int c = 0; c < n_components; ++c)
            {
              if (evaluation_flag & EvaluationFlags::values)
                eval0.template values<0, true, false>(values_dofs, values_quad);
              if (evaluation_flag & EvaluationFlags::gradients)
                eval0.template gradients<0, true, false>(values_dofs,
                                                         gradients_quad);

              // advance the next component in 1d array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += n_q_points;
            }
          break;

        case 2:
          for (unsigned int c = 0; c < n_components; ++c)
            {
              // grad x
              if (evaluation_flag & EvaluationFlags::gradients)
                {
                  eval0.template gradients<0, true, false>(values_dofs, temp1);
                  eval1.template values<1, true, false, 2>(temp1,
                                                           gradients_quad);
                }

              // grad y
              eval0.template values<0, true, false>(values_dofs, temp1);
              if (evaluation_flag & EvaluationFlags::gradients)
                eval1.template gradients<1, true, false, 2>(temp1,
                                                            gradients_quad + 1);

              // val: can use values applied in x
              if (evaluation_flag & EvaluationFlags::values)
                eval1.template values<1, true, false>(temp1, values_quad);

              // advance to the next component in 1d array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += 2 * n_q_points;
            }
          break;

        case 3:
          for (unsigned int c = 0; c < n_components; ++c)
            {
              if (evaluation_flag & EvaluationFlags::gradients)
                {
                  // grad x
                  eval0.template gradients<0, true, false>(values_dofs, temp1);
                  eval1.template values<1, true, false>(temp1, temp2);
                  eval2.template values<2, true, false, 3>(temp2,
                                                           gradients_quad);
                }

              // grad y
              eval0.template values<0, true, false>(values_dofs, temp1);
              if (evaluation_flag & EvaluationFlags::gradients)
                {
                  eval1.template gradients<1, true, false>(temp1, temp2);
                  eval2.template values<2, true, false, 3>(temp2,
                                                           gradients_quad + 1);
                }

              // grad z: can use the values applied in x direction stored in
              // temp1
              eval1.template values<1, true, false>(temp1, temp2);
              if (evaluation_flag & EvaluationFlags::gradients)
                eval2.template gradients<2, true, false, 3>(temp2,
                                                            gradients_quad + 2);

              // val: can use the values applied in x & y direction stored in
              // temp2
              if (evaluation_flag & EvaluationFlags::values)
                eval2.template values<2, true, false>(temp2, values_quad);

              // advance to the next component in 1d array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += 3 * n_q_points;
            }
          break;

        default:
          AssertThrow(false, ExcNotImplemented());
      }

    // case additional dof for FE_Q_DG0: add values; gradients and second
    // derivatives evaluate to zero
    if (type == MatrixFreeFunctions::tensor_symmetric_plus_dg0 &&
        (evaluation_flag & EvaluationFlags::values))
      {
        values_quad -= n_components * n_q_points;
        values_dofs -= n_components * dofs_per_comp;
        for (std::size_t c = 0; c < n_components; ++c)
          for (std::size_t q = 0; q < n_q_points; ++q)
            values_quad[c * n_q_points + q] +=
              values_dofs[(c + 1) * dofs_per_comp - 1];
      }
  }



  template <MatrixFreeFunctions::ElementType type,
            int                              dim,
            int                              fe_degree,
            int                              n_q_points_1d,
            typename Number>
  inline void
  FEEvaluationImpl<type, dim, fe_degree, n_q_points_1d, Number>::integrate(
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags integration_flag,
    Number                                *values_dofs_actual,
    FEEvaluationData<dim, Number, false>  &fe_eval,
    const bool                             add_into_values_array)
  {
    std::array<const MatrixFreeFunctions::UnivariateShapeData<Number2> *, 3>
      univariate_shape_data;

    const auto &shape_data = fe_eval.get_shape_info().data;
    univariate_shape_data.fill(&shape_data.front());

    if (shape_data.size() == dim)
      for (int i = 1; i < dim; ++i)
        univariate_shape_data[i] = &shape_data[i];

    Eval eval0 = create_evaluator_tensor_product(univariate_shape_data[0]);
    Eval eval1 = create_evaluator_tensor_product(univariate_shape_data[1]);
    Eval eval2 = create_evaluator_tensor_product(univariate_shape_data[2]);

    const unsigned int temp_size =
      Eval::n_rows_of_product == numbers::invalid_unsigned_int ?
        0 :
        (Eval::n_rows_of_product > Eval::n_columns_of_product ?
           Eval::n_rows_of_product :
           Eval::n_columns_of_product);
    Number *temp1 = fe_eval.get_scratch_data().begin();
    Number *temp2;
    if (temp_size == 0)
      {
        temp2 = temp1 + std::max(Utilities::fixed_power<dim>(
                                   shape_data.front().fe_degree + 1),
                                 Utilities::fixed_power<dim>(
                                   shape_data.front().n_q_points_1d));
      }
    else
      {
        temp2 = temp1 + temp_size;
      }

    const std::size_t  n_q_points = temp_size == 0 ?
                                      fe_eval.get_shape_info().n_q_points :
                                      Eval::n_columns_of_product;
    const unsigned int dofs_per_comp =
      (type == MatrixFreeFunctions::truncated_tensor) ?
        Utilities::fixed_power<dim>(shape_data.front().fe_degree + 1) :
        fe_eval.get_shape_info().dofs_per_component_on_cell;

    // expand dof_values to tensor product for truncated tensor products
    Number *values_dofs =
      (type == MatrixFreeFunctions::truncated_tensor) ?
        temp1 + 2 * (std::max<std::size_t>(
                      fe_eval.get_shape_info().dofs_per_component_on_cell,
                      n_q_points)) :
        values_dofs_actual;

    Number *values_quad    = fe_eval.begin_values();
    Number *gradients_quad = fe_eval.begin_gradients();

    switch (dim)
      {
        case 1:
          for (unsigned int c = 0; c < n_components; ++c)
            {
              if (integration_flag & EvaluationFlags::values)
                {
                  if (add_into_values_array == false)
                    eval0.template values<0, false, false>(values_quad,
                                                           values_dofs);
                  else
                    eval0.template values<0, false, true>(values_quad,
                                                          values_dofs);
                }
              if (integration_flag & EvaluationFlags::gradients)
                {
                  if (integration_flag & EvaluationFlags::values ||
                      add_into_values_array == true)
                    eval0.template gradients<0, false, true>(gradients_quad,
                                                             values_dofs);
                  else
                    eval0.template gradients<0, false, false>(gradients_quad,
                                                              values_dofs);
                }

              // advance to the next component in 1d array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += n_q_points;
            }
          break;

        case 2:
          for (unsigned int c = 0; c < n_components; ++c)
            {
              if ((integration_flag & EvaluationFlags::values) &&
                  !(integration_flag & EvaluationFlags::gradients))
                {
                  eval1.template values<1, false, false>(values_quad, temp1);
                  if (add_into_values_array == false)
                    eval0.template values<0, false, false>(temp1, values_dofs);
                  else
                    eval0.template values<0, false, true>(temp1, values_dofs);
                }
              if (integration_flag & EvaluationFlags::gradients)
                {
                  eval1.template gradients<1, false, false, 2>(gradients_quad +
                                                                 1,
                                                               temp1);
                  if (integration_flag & EvaluationFlags::values)
                    eval1.template values<1, false, true>(values_quad, temp1);
                  if (add_into_values_array == false)
                    eval0.template values<0, false, false>(temp1, values_dofs);
                  else
                    eval0.template values<0, false, true>(temp1, values_dofs);
                  eval1.template values<1, false, false, 2>(gradients_quad,
                                                            temp1);
                  eval0.template gradients<0, false, true>(temp1, values_dofs);
                }

              // advance to the next component in 1d array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += 2 * n_q_points;
            }
          break;

        case 3:
          for (unsigned int c = 0; c < n_components; ++c)
            {
              if ((integration_flag & EvaluationFlags::values) &&
                  !(integration_flag & EvaluationFlags::gradients))
                {
                  eval2.template values<2, false, false>(values_quad, temp1);
                  eval1.template values<1, false, false>(temp1, temp2);
                  if (add_into_values_array == false)
                    eval0.template values<0, false, false>(temp2, values_dofs);
                  else
                    eval0.template values<0, false, true>(temp2, values_dofs);
                }
              if (integration_flag & EvaluationFlags::gradients)
                {
                  eval2.template gradients<2, false, false, 3>(gradients_quad +
                                                                 2,
                                                               temp1);
                  if (integration_flag & EvaluationFlags::values)
                    eval2.template values<2, false, true>(values_quad, temp1);
                  eval1.template values<1, false, false>(temp1, temp2);
                  eval2.template values<2, false, false, 3>(gradients_quad + 1,
                                                            temp1);
                  eval1.template gradients<1, false, true>(temp1, temp2);
                  if (add_into_values_array == false)
                    eval0.template values<0, false, false>(temp2, values_dofs);
                  else
                    eval0.template values<0, false, true>(temp2, values_dofs);
                  eval2.template values<2, false, false, 3>(gradients_quad,
                                                            temp1);
                  eval1.template values<1, false, false>(temp1, temp2);
                  eval0.template gradients<0, false, true>(temp2, values_dofs);
                }

              // advance to the next component in 1d array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += 3 * n_q_points;
            }
          break;

        default:
          AssertThrow(false, ExcNotImplemented());
      }

    // case FE_Q_DG0: add values, gradients and second derivatives are zero
    if (type == MatrixFreeFunctions::tensor_symmetric_plus_dg0)
      {
        values_dofs -= n_components * dofs_per_comp - dofs_per_comp + 1;
        values_quad -= n_components * n_q_points;
        if (integration_flag & EvaluationFlags::values)
          for (unsigned int c = 0; c < n_components; ++c)
            {
              values_dofs[0] = values_quad[0];
              for (unsigned int q = 1; q < n_q_points; ++q)
                values_dofs[0] += values_quad[q];
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
            }
        else
          {
            for (unsigned int c = 0; c < n_components; ++c)
              values_dofs[c * dofs_per_comp] = Number();
            values_dofs += n_components * dofs_per_comp;
          }
      }

    if (type == MatrixFreeFunctions::truncated_tensor)
      truncate_tensor_product_to_complete_degrees<dim, fe_degree>(
        n_components,
        values_dofs_actual,
        values_dofs - dofs_per_comp * n_components,
        fe_eval);
  }



  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  inline void
  FEEvaluationImpl<
    MatrixFreeFunctions::tensor_none,
    dim,
    fe_degree,
    n_q_points_1d,
    Number>::evaluate(const unsigned int                     n_components,
                      const EvaluationFlags::EvaluationFlags evaluation_flag,
                      const Number                          *values_dofs_actual,
                      FEEvaluationData<dim, Number, false>  &fe_eval)
  {
    Assert(!(evaluation_flag & EvaluationFlags::hessians), ExcNotImplemented());

    const std::size_t n_dofs =
      fe_eval.get_shape_info().dofs_per_component_on_cell;
    const std::size_t n_q_points = fe_eval.get_shape_info().n_q_points;

    const auto &shape_data = fe_eval.get_shape_info().data;

    using Number2 =
      typename FEEvaluationData<dim, Number, false>::shape_info_number_type;

    if (evaluation_flag & EvaluationFlags::values)
      {
        const auto *const shape_values = shape_data.front().shape_values.data();
        auto             *out          = fe_eval.begin_values();
        const auto       *in           = values_dofs_actual;

        for (unsigned int c = 0; c < n_components; c += 3)
          {
            if (c + 1 == n_components)
              apply_matrix_vector_product<evaluate_general,
                                          EvaluatorQuantity::value,
                                          /*transpose_matrix*/ true,
                                          /*add*/ false,
                                          /*consider_strides*/ false,
                                          Number,
                                          Number2,
                                          /*n_components*/ 1>(
                shape_values, in, out, n_dofs, n_q_points, 1, 1);
            else if (c + 2 == n_components)
              apply_matrix_vector_product<evaluate_general,
                                          EvaluatorQuantity::value,
                                          /*transpose_matrix*/ true,
                                          /*add*/ false,
                                          /*consider_strides*/ false,
                                          Number,
                                          Number2,
                                          /*n_components*/ 2>(
                shape_values, in, out, n_dofs, n_q_points, 1, 1);
            else
              apply_matrix_vector_product<evaluate_general,
                                          EvaluatorQuantity::value,
                                          /*transpose_matrix*/ true,
                                          /*add*/ false,
                                          /*consider_strides*/ false,
                                          Number,
                                          Number2,
                                          /*n_components*/ 3>(
                shape_values, in, out, n_dofs, n_q_points, 1, 1);

            out += 3 * n_q_points;
            in += 3 * n_dofs;
          }
      }

    if (evaluation_flag & EvaluationFlags::gradients)
      {
        const auto *const shape_gradients =
          shape_data.front().shape_gradients.data();
        auto       *out = fe_eval.begin_gradients();
        const auto *in  = values_dofs_actual;

        for (unsigned int c = 0; c < n_components; c += 3)
          {
            if (c + 1 == n_components)
              apply_matrix_vector_product<evaluate_general,
                                          EvaluatorQuantity::value,
                                          /*transpose_matrix*/ true,
                                          /*add*/ false,
                                          /*consider_strides*/ false,
                                          Number,
                                          Number2,
                                          /*n_components*/ 1>(
                shape_gradients, in, out, n_dofs, n_q_points * dim, 1, 1);
            else if (c + 2 == n_components)
              apply_matrix_vector_product<evaluate_general,
                                          EvaluatorQuantity::value,
                                          /*transpose_matrix*/ true,
                                          /*add*/ false,
                                          /*consider_strides*/ false,
                                          Number,
                                          Number2,
                                          /*n_components*/ 2>(
                shape_gradients, in, out, n_dofs, n_q_points * dim, 1, 1);
            else
              apply_matrix_vector_product<evaluate_general,
                                          EvaluatorQuantity::value,
                                          /*transpose_matrix*/ true,
                                          /*add*/ false,
                                          /*consider_strides*/ false,
                                          Number,
                                          Number2,
                                          /*n_components*/ 3>(
                shape_gradients, in, out, n_dofs, n_q_points * dim, 1, 1);

            out += 3 * n_q_points * dim;
            in += 3 * n_dofs;
          }
      }
  }



  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  inline void
  FEEvaluationImpl<
    MatrixFreeFunctions::tensor_none,
    dim,
    fe_degree,
    n_q_points_1d,
    Number>::integrate(const unsigned int                     n_components,
                       const EvaluationFlags::EvaluationFlags integration_flag,
                       Number                               *values_dofs_actual,
                       FEEvaluationData<dim, Number, false> &fe_eval,
                       const bool add_into_values_array)
  {
    Assert(!(integration_flag & EvaluationFlags::hessians),
           ExcNotImplemented());

    const std::size_t n_dofs =
      fe_eval.get_shape_info().dofs_per_component_on_cell;
    const std::size_t n_q_points = fe_eval.get_shape_info().n_q_points;

    const auto &shape_data = fe_eval.get_shape_info().data;

    using Number2 =
      typename FEEvaluationData<dim, Number, false>::shape_info_number_type;

    if (integration_flag & EvaluationFlags::values)
      {
        const auto *const shape_values = shape_data.front().shape_values.data();
        auto             *in           = fe_eval.begin_values();
        auto             *out          = values_dofs_actual;

        for (unsigned int c = 0; c < n_components; c += 3)
          {
            if (add_into_values_array == false)
              {
                if (c + 1 == n_components)
                  apply_matrix_vector_product<evaluate_general,
                                              EvaluatorQuantity::value,
                                              /*transpose_matrix*/ false,
                                              /*add*/ false,
                                              /*consider_strides*/ false,
                                              Number,
                                              Number2,
                                              /*n_components*/ 1>(
                    shape_values, in, out, n_dofs, n_q_points, 1, 1);
                else if (c + 2 == n_components)
                  apply_matrix_vector_product<evaluate_general,
                                              EvaluatorQuantity::value,
                                              /*transpose_matrix*/ false,
                                              /*add*/ false,
                                              /*consider_strides*/ false,
                                              Number,
                                              Number2,
                                              /*n_components*/ 2>(
                    shape_values, in, out, n_dofs, n_q_points, 1, 1);
                else
                  apply_matrix_vector_product<evaluate_general,
                                              EvaluatorQuantity::value,
                                              /*transpose_matrix*/ false,
                                              /*add*/ false,
                                              /*consider_strides*/ false,
                                              Number,
                                              Number2,
                                              /*n_components*/ 3>(
                    shape_values, in, out, n_dofs, n_q_points, 1, 1);
              }
            else
              {
                if (c + 1 == n_components)
                  apply_matrix_vector_product<evaluate_general,
                                              EvaluatorQuantity::value,
                                              /*transpose_matrix*/ false,
                                              /*add*/ true,
                                              /*consider_strides*/ false,
                                              Number,
                                              Number2,
                                              /*n_components*/ 1>(
                    shape_values, in, out, n_dofs, n_q_points, 1, 1);
                else if (c + 2 == n_components)
                  apply_matrix_vector_product<evaluate_general,
                                              EvaluatorQuantity::value,
                                              /*transpose_matrix*/ false,
                                              /*add*/ true,
                                              /*consider_strides*/ false,
                                              Number,
                                              Number2,
                                              /*n_components*/ 2>(
                    shape_values, in, out, n_dofs, n_q_points, 1, 1);
                else
                  apply_matrix_vector_product<evaluate_general,
                                              EvaluatorQuantity::value,
                                              /*transpose_matrix*/ false,
                                              /*add*/ true,
                                              /*consider_strides*/ false,
                                              Number,
                                              Number2,
                                              /*n_components*/ 3>(
                    shape_values, in, out, n_dofs, n_q_points, 1, 1);
              }
            out += 3 * n_dofs;
            in += 3 * n_q_points;
          }
      }

    if (integration_flag & EvaluationFlags::gradients)
      {
        const auto *const shape_gradients =
          shape_data.front().shape_gradients.data();
        auto *in  = fe_eval.begin_gradients();
        auto *out = values_dofs_actual;

        for (unsigned int c = 0; c < n_components; c += 3)
          {
            if (add_into_values_array == false &&
                !(integration_flag & EvaluationFlags::values))
              {
                if (c + 1 == n_components)
                  apply_matrix_vector_product<evaluate_general,
                                              EvaluatorQuantity::value,
                                              /*transpose_matrix*/ false,
                                              /*add*/ false,
                                              /*consider_strides*/ false,
                                              Number,
                                              Number2,
                                              /*n_components*/ 1>(
                    shape_gradients, in, out, n_dofs, n_q_points * dim, 1, 1);
                else if (c + 2 == n_components)
                  apply_matrix_vector_product<evaluate_general,
                                              EvaluatorQuantity::value,
                                              /*transpose_matrix*/ false,
                                              /*add*/ false,
                                              /*consider_strides*/ false,
                                              Number,
                                              Number2,
                                              /*n_components*/ 2>(
                    shape_gradients, in, out, n_dofs, n_q_points * dim, 1, 1);
                else
                  apply_matrix_vector_product<evaluate_general,
                                              EvaluatorQuantity::value,
                                              /*transpose_matrix*/ false,
                                              /*add*/ false,
                                              /*consider_strides*/ false,
                                              Number,
                                              Number2,
                                              /*n_components*/ 3>(
                    shape_gradients, in, out, n_dofs, n_q_points * dim, 1, 1);
              }
            else
              {
                if (c + 1 == n_components)
                  apply_matrix_vector_product<evaluate_general,
                                              EvaluatorQuantity::value,
                                              /*transpose_matrix*/ false,
                                              /*add*/ true,
                                              /*consider_strides*/ false,
                                              Number,
                                              Number2,
                                              /*n_components*/ 1>(
                    shape_gradients, in, out, n_dofs, n_q_points * dim, 1, 1);
                else if (c + 2 == n_components)
                  apply_matrix_vector_product<evaluate_general,
                                              EvaluatorQuantity::value,
                                              /*transpose_matrix*/ false,
                                              /*add*/ true,
                                              /*consider_strides*/ false,
                                              Number,
                                              Number2,
                                              /*n_components*/ 2>(
                    shape_gradients, in, out, n_dofs, n_q_points * dim, 1, 1);
                else
                  apply_matrix_vector_product<evaluate_general,
                                              EvaluatorQuantity::value,
                                              /*transpose_matrix*/ false,
                                              /*add*/ true,
                                              /*consider_strides*/ false,
                                              Number,
                                              Number2,
                                              /*n_components*/ 3>(
                    shape_gradients, in, out, n_dofs, n_q_points * dim, 1, 1);
              }
            out += 3 * n_dofs;
            in += 3 * n_q_points * dim;
          }
      }
  }



  /**
   * This struct implements the change between two different bases. This is an
   * ingredient in the FEEvaluationImplTransformToCollocation class where we
   * first transform to the appropriate basis where we can compute the
   * derivative through collocation techniques.
   *
   * This class allows for dimension-independent application of the operation,
   * implemented by template recursion. It has been tested up to 6d.
   */
  template <EvaluatorVariant  variant,
            EvaluatorQuantity quantity,
            int               dim,
            int               basis_size_1,
            int               basis_size_2>
  struct FEEvaluationImplBasisChange
  {
    static_assert(basis_size_1 == 0 || basis_size_1 <= basis_size_2,
                  "The second dimension must not be smaller than the first");

    /**
     * This applies the transformation that contracts over the rows of the
     * coefficient array, generating values along the columns of the
     * coefficient array.
     *
     * @param n_components The number of vector components.
     * @param transformation_matrix The coefficient matrix handed in as a
     *                     vector, using @p basis_size_1 rows and @p basis_size_2
     *                     columns if interpreted as a matrix.
     * @param values_in    The array of the input of size basis_size_1^dim. It
     *                     may alias with values_out
     * @param values_out   The array of size basis_size_2^dim where the results
     *                     of the transformation are stored. It may alias with
     *                     the values_in array.
     * @param basis_size_1_variable In case the template argument
     * @p basis_size_1 is zero, the size of the first basis can alternatively
     * be passed in as a run time argument. The template argument takes
     * precedence in case it is nonzero for efficiency reasons.
     * @param basis_size_2_variable In case the template argument
     * @p basis_size_1 is zero, the size of the second basis can alternatively
     * be passed in as a run time argument.
     */
    template <typename Number, typename Number2>
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
      static void
      do_forward(const unsigned int            n_components,
                 const AlignedVector<Number2> &transformation_matrix,
                 const Number                 *values_in,
                 Number                       *values_out,
                 const unsigned int            basis_size_1_variable =
                   numbers::invalid_unsigned_int,
                 const unsigned int basis_size_2_variable =
                   numbers::invalid_unsigned_int)
    {
      Assert(
        basis_size_1 != 0 || basis_size_1_variable <= basis_size_2_variable,
        ExcMessage("The second dimension must not be smaller than the first"));

      Assert(quantity == EvaluatorQuantity::value, ExcInternalError());

      // we do recursion until dim==1 or dim==2 and we have
      // basis_size_1==basis_size_2. The latter optimization increases
      // optimization possibilities for the compiler but does only work for
      // aliased pointers if the sizes are equal.
      constexpr int next_dim = (dim == 1 || (dim == 2 && basis_size_1 > 0 &&
                                             basis_size_1 == basis_size_2)) ?
                                 dim :
                                 dim - 1;

      EvaluatorTensorProduct<variant,
                             dim,
                             basis_size_1,
                             (basis_size_1 == 0 ? 0 : basis_size_2),
                             Number,
                             Number2>
                         eval_val(transformation_matrix,
                                  {},
                                  {},
                 basis_size_1_variable,
                 basis_size_2_variable);
      const unsigned int np_1 =
        basis_size_1 > 0 ? basis_size_1 : basis_size_1_variable;
      const unsigned int np_2 =
        basis_size_1 > 0 ? basis_size_2 : basis_size_2_variable;
      Assert(np_1 > 0 && np_1 != numbers::invalid_unsigned_int,
             ExcMessage("Cannot transform with 0-point basis"));
      Assert(np_2 > 0 && np_2 != numbers::invalid_unsigned_int,
             ExcMessage("Cannot transform with 0-point basis"));

      // run loop backwards to ensure correctness if values_in aliases with
      // values_out in case with basis_size_1 < basis_size_2
      values_in = values_in + n_components * Utilities::fixed_power<dim>(np_1);
      values_out =
        values_out + n_components * Utilities::fixed_power<dim>(np_2);
      for (unsigned int c = n_components; c != 0; --c)
        {
          values_in -= Utilities::fixed_power<dim>(np_1);
          values_out -= Utilities::fixed_power<dim>(np_2);
          if (next_dim < dim)
            for (unsigned int q = np_1; q != 0; --q)
              FEEvaluationImplBasisChange<variant,
                                          quantity,
                                          next_dim,
                                          basis_size_1,
                                          basis_size_2>::
                do_forward(1,
                           transformation_matrix,
                           values_in +
                             (q - 1) * Utilities::fixed_power<next_dim>(np_1),
                           values_out +
                             (q - 1) * Utilities::fixed_power<next_dim>(np_2),
                           basis_size_1_variable,
                           basis_size_2_variable);

          // the recursion stops if dim==1 or if dim==2 and
          // basis_size_1==basis_size_2 (the latter is used because the
          // compiler generates nicer code)
          if (basis_size_1 > 0 && basis_size_2 == basis_size_1 && dim == 2)
            {
              eval_val.template values<0, true, false>(values_in, values_out);
              eval_val.template values<1, true, false>(values_out, values_out);
            }
          else if (dim == 1)
            eval_val.template values<dim - 1, true, false>(values_in,
                                                           values_out);
          else
            eval_val.template values<dim - 1, true, false>(values_out,
                                                           values_out);
        }
    }

    /**
     * This applies the transformation that contracts over the columns of the
     * coefficient array, generating values along the rows of the coefficient
     * array.
     *
     * @param n_components The number of vector components.
     * @param transformation_matrix The coefficient matrix handed in as a
     *                     vector, using @p basis_size_1 rows and @p basis_size_2
     *                     columns if interpreted as a matrix.
     * @param add_into_result Define whether the result should be added into the
     *                     array @p values_out (if true) or overwrite the
     *                     previous content. The result is undefined in case
     *                     values_in and values_out point to the same array and
     *                     @p add_into_result is true, in which case an
     *                     exception is thrown.
     * @param values_in    The array of the input of size basis_size_2^dim. It
     *                     may alias with values_out. Note that the previous
     *                     content of @p values_in is overwritten within the
     *                     function.
     * @param values_out   The array of size basis_size_1^dim where the results
     *                     of the transformation are stored. It may alias with
     *                     the @p values_in array.
     * @param basis_size_1_variable In case the template argument
     * @p basis_size_1 is zero, the size of the first basis can alternatively
     * be passed in as a run time argument. The template argument takes
     * precedence in case it is nonzero for efficiency reasons.
     * @param basis_size_2_variable In case the template argument
     * @p basis_size_1 is zero, the size of the second basis can alternatively
     * be passed in as a run time argument.
     */
    template <typename Number, typename Number2>
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
      static void
      do_backward(const unsigned int            n_components,
                  const AlignedVector<Number2> &transformation_matrix,
                  const bool                    add_into_result,
                  Number                       *values_in,
                  Number                       *values_out,
                  const unsigned int            basis_size_1_variable =
                    numbers::invalid_unsigned_int,
                  const unsigned int basis_size_2_variable =
                    numbers::invalid_unsigned_int)
    {
      Assert(
        basis_size_1 != 0 || basis_size_1_variable <= basis_size_2_variable,
        ExcMessage("The second dimension must not be smaller than the first"));
      Assert(add_into_result == false || values_in != values_out,
             ExcMessage(
               "Input and output cannot alias with each other when "
               "adding the result of the basis change to existing data"));

      Assert(quantity == EvaluatorQuantity::value ||
               quantity == EvaluatorQuantity::hessian,
             ExcInternalError());

      constexpr int next_dim =
        (dim > 2 ||
         ((basis_size_1 == 0 || basis_size_2 > basis_size_1) && dim > 1)) ?
          dim - 1 :
          dim;
      EvaluatorTensorProduct<variant,
                             dim,
                             basis_size_1,
                             (basis_size_1 == 0 ? 0 : basis_size_2),
                             Number,
                             Number2>
                         eval_val(transformation_matrix,
                 transformation_matrix,
                 transformation_matrix,
                 basis_size_1_variable,
                 basis_size_2_variable);
      const unsigned int np_1 =
        basis_size_1 > 0 ? basis_size_1 : basis_size_1_variable;
      const unsigned int np_2 =
        basis_size_1 > 0 ? basis_size_2 : basis_size_2_variable;
      Assert(np_1 > 0 && np_1 != numbers::invalid_unsigned_int,
             ExcMessage("Cannot transform with 0-point basis"));
      Assert(np_2 > 0 && np_2 != numbers::invalid_unsigned_int,
             ExcMessage("Cannot transform with 0-point basis"));

      for (unsigned int c = 0; c < n_components; ++c)
        {
          if (basis_size_1 > 0 && basis_size_2 == basis_size_1 && dim == 2)
            {
              if (quantity == EvaluatorQuantity::value)
                eval_val.template values<1, false, false>(values_in, values_in);
              else
                eval_val.template hessians<1, false, false>(values_in,
                                                            values_in);

              if (add_into_result)
                {
                  if (quantity == EvaluatorQuantity::value)
                    eval_val.template values<0, false, true>(values_in,
                                                             values_out);
                  else
                    eval_val.template hessians<0, false, true>(values_in,
                                                               values_out);
                }
              else
                {
                  if (quantity == EvaluatorQuantity::value)
                    eval_val.template values<0, false, false>(values_in,
                                                              values_out);
                  else
                    eval_val.template hessians<0, false, false>(values_in,
                                                                values_out);
                }
            }
          else
            {
              if (dim == 1 && add_into_result)
                {
                  if (quantity == EvaluatorQuantity::value)
                    eval_val.template values<0, false, true>(values_in,
                                                             values_out);
                  else
                    eval_val.template hessians<0, false, true>(values_in,
                                                               values_out);
                }
              else if (dim == 1)
                {
                  if (quantity == EvaluatorQuantity::value)
                    eval_val.template values<0, false, false>(values_in,
                                                              values_out);
                  else
                    eval_val.template hessians<0, false, false>(values_in,
                                                                values_out);
                }
              else
                {
                  if (quantity == EvaluatorQuantity::value)
                    eval_val.template values<dim - 1, false, false>(values_in,
                                                                    values_in);
                  else
                    eval_val.template hessians<dim - 1, false, false>(
                      values_in, values_in);
                }
            }
          if (next_dim < dim)
            for (unsigned int q = 0; q < np_1; ++q)
              FEEvaluationImplBasisChange<variant,
                                          quantity,
                                          next_dim,
                                          basis_size_1,
                                          basis_size_2>::
                do_backward(1,
                            transformation_matrix,
                            add_into_result,
                            values_in +
                              q * Utilities::fixed_power<next_dim>(np_2),
                            values_out +
                              q * Utilities::fixed_power<next_dim>(np_1),
                            basis_size_1_variable,
                            basis_size_2_variable);

          values_in += Utilities::fixed_power<dim>(np_2);
          values_out += Utilities::fixed_power<dim>(np_1);
        }
    }

    /**
     * This operation applies a mass-matrix-like operation, consisting of a
     * do_forward() operation, multiplication by the coefficients in the
     * quadrature points, and the do_backward() operation.
     *
     * @param n_components The number of vector components.
     * @param transformation_matrix The coefficient matrix handed in as a
     *                     vector, using @p basis_size_1 rows and @p basis_size_2
     *                     columns if interpreted as a matrix.
     * @param coefficients The array of coefficients by which the result is
     *                     multiplied. Its length must be either
     *                     basis_size_2^dim or n_components*basis_size_2^dim.
     * @param values_in    The array of the input of size basis_size_2^dim. It
     *                     may alias with values_out.
     * @param scratch_data Array to hold temporary data during the operation.
     *                     Must be of length basis_size_2^dim.
     * @param values_out   The array of size basis_size_1^dim where the results
     *                     of the transformation are stored. It may alias with
     *                     the values_in array.
     */
    template <typename Number, typename Number2>
    static void
    do_mass(const unsigned int            n_components,
            const AlignedVector<Number2> &transformation_matrix,
            const AlignedVector<Number>  &coefficients,
            const Number                 *values_in,
            Number                       *scratch_data,
            Number                       *values_out)
    {
      constexpr int next_dim = dim > 1 ? dim - 1 : dim;
      Number       *my_scratch =
        basis_size_1 != basis_size_2 ? scratch_data : values_out;

      const unsigned int size_per_component = Utilities::pow(basis_size_2, dim);
      Assert(coefficients.size() == size_per_component ||
               coefficients.size() == n_components * size_per_component,
             ExcDimensionMismatch(coefficients.size(), size_per_component));
      const unsigned int stride =
        coefficients.size() == size_per_component ? 0 : 1;

      for (unsigned int q = basis_size_1; q != 0; --q)
        FEEvaluationImplBasisChange<
          variant,
          EvaluatorQuantity::value,
          next_dim,
          basis_size_1,
          basis_size_2>::do_forward(n_components,
                                    transformation_matrix,
                                    values_in +
                                      (q - 1) *
                                        Utilities::pow(basis_size_1, dim - 1),
                                    my_scratch +
                                      (q - 1) *
                                        Utilities::pow(basis_size_2, dim - 1));
      EvaluatorTensorProduct<variant,
                             dim,
                             basis_size_1,
                             basis_size_2,
                             Number,
                             Number2>
                         eval_val(transformation_matrix);
      const unsigned int n_inner_blocks =
        (dim > 1 && basis_size_2 < 10) ? basis_size_2 : 1;
      const unsigned int n_blocks = Utilities::pow(basis_size_2, dim - 1);
      for (unsigned int ii = 0; ii < n_blocks; ii += n_inner_blocks)
        for (unsigned int c = 0; c < n_components; ++c)
          {
            for (unsigned int i = ii; i < ii + n_inner_blocks; ++i)
              eval_val.template values_one_line<dim - 1, true, false>(
                my_scratch + i, my_scratch + i);
            for (unsigned int q = 0; q < basis_size_2; ++q)
              for (unsigned int i = ii; i < ii + n_inner_blocks; ++i)
                my_scratch[i + q * n_blocks + c * size_per_component] *=
                  coefficients[i + q * n_blocks +
                               c * stride * size_per_component];
            for (unsigned int i = ii; i < ii + n_inner_blocks; ++i)
              eval_val.template values_one_line<dim - 1, false, false>(
                my_scratch + i, my_scratch + i);
          }
      for (unsigned int q = 0; q < basis_size_1; ++q)
        FEEvaluationImplBasisChange<variant,
                                    EvaluatorQuantity::value,
                                    next_dim,
                                    basis_size_1,
                                    basis_size_2>::
          do_backward(n_components,
                      transformation_matrix,
                      false,
                      my_scratch + q * Utilities::pow(basis_size_2, dim - 1),
                      values_out + q * Utilities::pow(basis_size_1, dim - 1));
    }
  };



  /**
   * Internal function that evaluates the gradients of finite element
   * functions represented by bases in the collocation space, used by
   * FEEvaluationImplCollocation and FEEvaluationImplTransformToCollocation.
   * The evaluation strategy uses sum factorization with the even-odd
   * optimization and fixed loop bounds.
   */
  template <int n_points_1d, int dim, typename Number, typename Number2>
  inline void
  evaluate_gradients_collocation(
    const MatrixFreeFunctions::UnivariateShapeData<Number2> &shape,
    const Number                                            *values,
    Number                                                  *gradients)
  {
    AssertDimension(shape.shape_gradients_collocation_eo.size(),
                    (n_points_1d + 1) / 2 * n_points_1d);

    EvaluatorTensorProduct<evaluate_evenodd,
                           dim,
                           n_points_1d,
                           n_points_1d,
                           Number,
                           Number2>
      eval({}, shape.shape_gradients_collocation_eo, {});
    EvaluatorTensorProduct<evaluate_evenodd,
                           2,
                           n_points_1d,
                           n_points_1d,
                           Number,
                           Number2>
      eval_2d({}, shape.shape_gradients_collocation_eo, {});

    if (dim == 1)
      eval.template gradients<0, true, false>(values, gradients);
    else
      {
        if (dim > 2)
          eval.template gradients<2, true, false, dim>(values, gradients + 2);
        constexpr unsigned int loop_bound  = (dim > 2 ? n_points_1d : 1);
        constexpr unsigned int n_points_2d = n_points_1d * n_points_1d;
        const Number          *in = values + (loop_bound - 1) * n_points_2d;
        Number *out = gradients + (loop_bound - 1) * dim * n_points_2d;
        for (unsigned int l = 0; l < loop_bound; ++l)
          {
            eval_2d.template gradients<0, true, false, dim>(in, out);
            eval_2d.template gradients<1, true, false, dim>(in, out + 1);
            in -= n_points_2d;
            out -= dim * n_points_2d;
          }
      }
  }



  /**
   * Internal function that multiplies by the gradients of test functions and
   * sums over quadrature points for function representations in the
   * collocation space, used by FEEvaluationImplCollocation and
   * FEEvaluationImplTransformToCollocation. The evaluation strategy uses sum
   * factorization with the even-odd optimization and fixed loop bounds.
   */
  template <int n_points_1d, int dim, typename Number, typename Number2>
  inline void
  integrate_gradients_collocation(
    const MatrixFreeFunctions::UnivariateShapeData<Number2> &shape,
    Number                                                  *values,
    const Number                                            *gradients,
    const bool add_into_values_array)
  {
    AssertDimension(shape.shape_gradients_collocation_eo.size(),
                    (n_points_1d + 1) / 2 * n_points_1d);

    EvaluatorTensorProduct<evaluate_evenodd,
                           dim,
                           n_points_1d,
                           n_points_1d,
                           Number,
                           Number2>
      eval({}, shape.shape_gradients_collocation_eo, {});
    EvaluatorTensorProduct<evaluate_evenodd,
                           2,
                           n_points_1d,
                           n_points_1d,
                           Number,
                           Number2>
      eval_2d({}, shape.shape_gradients_collocation_eo, {});

    if (dim == 1)
      {
        if (add_into_values_array)
          eval.template gradients<0, false, true>(gradients, values);
        else
          eval.template gradients<0, false, false>(gradients, values);
      }
    else
      {
        constexpr unsigned int loop_bound  = (dim > 2 ? n_points_1d : 1);
        constexpr unsigned int n_points_2d = n_points_1d * n_points_1d;

        const Number *in  = gradients + (loop_bound - 1) * dim * n_points_2d;
        Number       *out = values + (loop_bound - 1) * n_points_2d;
        for (unsigned int l = 0; l < loop_bound; ++l)
          {
            if (add_into_values_array)
              eval_2d.template gradients<0, false, true, dim>(in, out);
            else
              eval_2d.template gradients<0, false, false, dim>(in, out);
            eval_2d.template gradients<1, false, true, dim>(in + 1, out);
            in -= dim * n_points_2d;
            out -= n_points_2d;
          }
      }
    if (dim > 2)
      eval.template gradients<2, false, true, dim>(gradients + 2, values);
  }



  /**
   * Internal function that evaluates the Hessians of finite element functions
   * represented by bases in the collocation space, used by
   * FEEvaluationImplSelector. The evaluation strategy uses sum
   * factorization with fixed loop bounds.
   */
  template <int n_points_1d, int dim, typename Number>
  inline void
  evaluate_hessians_collocation(const unsigned int n_components,
                                FEEvaluationData<dim, Number, false> &fe_eval)
  {
    using Number2 =
      typename FEEvaluationData<dim, Number, false>::shape_info_number_type;

    // might have non-symmetric quadrature formula, so use the more
    // conservative 'evaluate_general' scheme rather than 'even_odd' as the
    // Hessians are not used very often
    const MatrixFreeFunctions::UnivariateShapeData<Number2> &data =
      fe_eval.get_shape_info().data[0];
    AssertDimension(data.shape_gradients_collocation.size(),
                    data.n_q_points_1d * data.n_q_points_1d);
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           n_points_1d,
                           n_points_1d,
                           Number,
                           Number2>
      eval({},
           data.shape_gradients_collocation.data(),
           data.shape_hessians_collocation.data(),
           data.n_q_points_1d,
           data.n_q_points_1d);

    const Number     *values   = fe_eval.begin_values();
    Number           *hessians = fe_eval.begin_hessians();
    Number           *scratch  = fe_eval.get_scratch_data().begin();
    const std::size_t n_points = fe_eval.get_shape_info().n_q_points;
    for (unsigned int comp = 0; comp < n_components; ++comp)
      {
        // xx derivative
        eval.template hessians<0, true, false>(values, hessians);
        if (dim > 1)
          {
            // xy derivative: we might or might not have the gradients already
            // computed elsewhere, but we recompute them here since it adds
            // only moderate extra work (at most 25%)
            eval.template gradients<0, true, false>(values, scratch);
            eval.template gradients<1, true, false>(scratch,
                                                    hessians + dim * n_points);
            // yy derivative
            eval.template hessians<1, true, false>(values, hessians + n_points);
          }
        if (dim > 2)
          {
            // xz derivative
            eval.template gradients<2, true, false>(scratch,
                                                    hessians + 4 * n_points);
            // yz derivative
            eval.template gradients<1, true, false>(values, scratch);
            eval.template gradients<2, true, false>(scratch,
                                                    hessians + 5 * n_points);
            // zz derivative
            eval.template hessians<2, true, false>(values,
                                                   hessians + 2 * n_points);
          }

        values += n_points;
        hessians += (dim * (dim + 1)) / 2 * n_points;
      }
  }



  /**
   * Internal function that multiplies by the Hessians of test functions and
   * sums over quadrature points for function representations in the
   * collocation space, used by FEEvaluationImplSelector. The evaluation
   * strategy uses sum factorization with fixed loop bounds.
   */
  template <int n_q_points_1d, int dim, typename Number>
  inline void
  integrate_hessians_collocation(const unsigned int n_components,
                                 FEEvaluationData<dim, Number, false> &fe_eval,
                                 const bool add_into_values_array)
  {
    using Number2 =
      typename FEEvaluationData<dim, Number, false>::shape_info_number_type;

    const MatrixFreeFunctions::UnivariateShapeData<Number2> &data =
      fe_eval.get_shape_info().data[0];
    AssertDimension(data.shape_gradients_collocation.size(),
                    data.n_q_points_1d * data.n_q_points_1d);
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           n_q_points_1d,
                           n_q_points_1d,
                           Number,
                           Number2>
                      eval({},
           data.shape_gradients_collocation.data(),
           data.shape_hessians_collocation.data(),
           data.n_q_points_1d,
           data.n_q_points_1d);
    Number           *values   = fe_eval.begin_values();
    const Number     *hessians = fe_eval.begin_hessians();
    Number           *scratch  = fe_eval.get_scratch_data().begin();
    const std::size_t n_points = fe_eval.get_shape_info().n_q_points;

    for (unsigned int comp = 0; comp < n_components; ++comp)
      {
        // xx derivative
        if (add_into_values_array == true)
          eval.template hessians<0, false, true>(hessians, values);
        else
          eval.template hessians<0, false, false>(hessians, values);

        // yy derivative
        if (dim > 1)
          eval.template hessians<1, false, true>(hessians + n_points, values);
        if (dim > 2)
          {
            // zz derivative
            eval.template hessians<2, false, true>(hessians + 2 * n_points,
                                                   values);
            // yz derivative
            eval.template gradients<2, false, false>(hessians + 5 * n_points,
                                                     scratch);
            eval.template gradients<1, false, true>(scratch, values);

            // xz derivative
            eval.template gradients<2, false, false>(hessians + 4 * n_points,
                                                     scratch);
          }

        if (dim > 1)
          {
            // xy derivative, combined with xz in 3d
            eval.template gradients<1, false, (dim > 2)>(hessians +
                                                           dim * n_points,
                                                         scratch);
            eval.template gradients<0, false, true>(scratch, values);
          }

        values += n_points;
        hessians += (dim * (dim + 1)) / 2 * n_points;
      }
  }



  /**
   * Internal function to evaluate the Hessians of finite element functions in
   * the non-collocation setting as a fall-back. The evaluation strategy uses
   * sum factorization with run-time loop bounds and is thus slower than the
   * collocation case, but it is not as widely used and thus uncritical.
   */
  template <int dim, typename Number>
  void
  evaluate_hessians_slow(const unsigned int                    n_components,
                         const Number                         *values_dofs,
                         FEEvaluationData<dim, Number, false> &fe_eval)
  {
    const auto &univariate_shape_data = fe_eval.get_shape_info().data;
    using Impl =
      FEEvaluationImpl<MatrixFreeFunctions::tensor_general, dim, -1, 0, Number>;
    using Eval = typename Impl::Eval;
    Eval eval0 =
      Impl::create_evaluator_tensor_product(&univariate_shape_data[0]);
    Eval eval1 = Impl::create_evaluator_tensor_product(
      &univariate_shape_data[std::min<int>(1,
                                           univariate_shape_data.size() - 1)]);
    Eval eval2 = Impl::create_evaluator_tensor_product(
      &univariate_shape_data[std::min<int>(2,
                                           univariate_shape_data.size() - 1)]);

    const unsigned int n_points = fe_eval.get_shape_info().n_q_points;
    Number            *tmp1     = fe_eval.get_scratch_data().begin();
    Number            *tmp2 =
      tmp1 + std::max(Utilities::fixed_power<dim>(
                        univariate_shape_data.front().fe_degree + 1),
                      Utilities::fixed_power<dim>(
                        univariate_shape_data.front().n_q_points_1d));
    Number *hessians = fe_eval.begin_hessians();

    for (unsigned int comp = 0; comp < n_components;
         ++comp,
                      hessians += n_points * dim * (dim + 1) / 2,
                      values_dofs +=
                      fe_eval.get_shape_info().dofs_per_component_on_cell)
      switch (dim)
        {
          case 1:
            eval0.template hessians<0, true, false>(values_dofs, hessians);
            break;
          case 2:
            // xx derivative
            eval0.template hessians<0, true, false>(values_dofs, tmp1);
            eval1.template values<1, true, false>(tmp1, hessians);
            // xy derivative
            eval0.template gradients<0, true, false>(values_dofs, tmp1);
            eval1.template gradients<1, true, false>(tmp1,
                                                     hessians + 2 * n_points);
            // yy derivative
            eval0.template values<0, true, false>(values_dofs, tmp1);
            eval1.template hessians<1, true, false>(tmp1, hessians + n_points);
            break;
          case 3:
            // xx derivative
            eval0.template hessians<0, true, false>(values_dofs, tmp1);
            eval1.template values<1, true, false>(tmp1, tmp2);
            eval2.template values<2, true, false>(tmp2, hessians);
            // xy derivative
            eval0.template gradients<0, true, false>(values_dofs, tmp1);
            eval1.template gradients<1, true, false>(tmp1, tmp2);
            eval2.template values<2, true, false>(tmp2,
                                                  hessians + 3 * n_points);
            // xz derivative
            eval1.template values<1, true, false>(tmp1, tmp2);
            eval2.template gradients<2, true, false>(tmp2,
                                                     hessians + 4 * n_points);
            // yy derivative
            eval0.template values<0, true, false>(values_dofs, tmp1);
            eval1.template hessians<1, true, false>(tmp1, tmp2);
            eval2.template values<2, true, false>(tmp2, hessians + n_points);
            // yz derivative
            eval1.template gradients<1, true, false>(tmp1, tmp2);
            eval2.template gradients<2, true, false>(tmp2,
                                                     hessians + 5 * n_points);
            // zz derivative
            eval1.template values<1, true, false>(tmp1, tmp2);
            eval2.template hessians<2, true, false>(tmp2,
                                                    hessians + 2 * n_points);
            break;

          default:
            Assert(false,
                   ExcNotImplemented(
                     "Only 1d, 2d and 3d implemented for Hessian"));
        }
  }



  /**
   * Internal function to multiply by the Hessians of the test functions and
   * integrate in the non-collocation setting as a fall-back. The evaluation
   * strategy uses sum factorization with run-time loop bounds and is thus
   * slower than the collocation case, but it is not as widely used and thus
   * uncritical.
   */
  template <int dim, typename Number>
  void
  integrate_hessians_slow(const unsigned int n_components,
                          const FEEvaluationData<dim, Number, false> &fe_eval,
                          Number    *values_dofs,
                          const bool add_into_values_array)
  {
    const auto &univariate_shape_data = fe_eval.get_shape_info().data;
    using Impl =
      FEEvaluationImpl<MatrixFreeFunctions::tensor_general, dim, -1, 0, Number>;
    using Eval = typename Impl::Eval;
    Eval eval0 =
      Impl::create_evaluator_tensor_product(&univariate_shape_data[0]);
    Eval eval1 = Impl::create_evaluator_tensor_product(
      &univariate_shape_data[std::min<int>(1,
                                           univariate_shape_data.size() - 1)]);
    Eval eval2 = Impl::create_evaluator_tensor_product(
      &univariate_shape_data[std::min<int>(2,
                                           univariate_shape_data.size() - 1)]);

    const unsigned int n_points = fe_eval.get_shape_info().n_q_points;
    Number            *tmp1     = fe_eval.get_scratch_data().begin();
    Number            *tmp2 =
      tmp1 + std::max(Utilities::fixed_power<dim>(
                        univariate_shape_data.front().fe_degree + 1),
                      Utilities::fixed_power<dim>(
                        univariate_shape_data.front().n_q_points_1d));
    const Number *hessians = fe_eval.begin_hessians();

    for (unsigned int comp = 0; comp < n_components;
         ++comp,
                      hessians += n_points * dim * (dim + 1) / 2,
                      values_dofs +=
                      fe_eval.get_shape_info().dofs_per_component_on_cell)
      switch (dim)
        {
          case 1:
            if (add_into_values_array)
              eval0.template hessians<0, false, true>(hessians, values_dofs);
            else
              eval0.template hessians<0, false, false>(hessians, values_dofs);
            break;
          case 2:
            // xx derivative
            eval1.template values<1, false, false>(hessians, tmp1);
            if (add_into_values_array)
              eval0.template hessians<0, false, true>(tmp1, values_dofs);
            else
              eval0.template hessians<0, false, false>(tmp1, values_dofs);

            // xy derivative
            eval1.template gradients<1, false, false>(hessians + 2 * n_points,
                                                      tmp1);
            eval0.template gradients<0, false, true>(tmp1, values_dofs);
            // yy derivative
            eval1.template hessians<1, false, false>(hessians + n_points, tmp1);
            eval0.template values<0, false, true>(tmp1, values_dofs);
            break;
          case 3:
            // xx derivative
            eval2.template values<2, false, false>(hessians, tmp1);
            eval1.template values<1, false, false>(tmp1, tmp2);

            if (add_into_values_array)
              eval0.template hessians<0, false, true>(tmp2, values_dofs);
            else
              eval0.template hessians<0, false, false>(tmp2, values_dofs);

            // xy derivative
            eval2.template values<2, false, false>(hessians + 3 * n_points,
                                                   tmp1);
            eval1.template gradients<1, false, false>(tmp1, tmp2);
            // xz derivative
            eval2.template gradients<2, false, false>(hessians + 4 * n_points,
                                                      tmp1);
            eval1.template values<1, false, true>(tmp1, tmp2);
            eval1.template values<0, false, true>(tmp2, values_dofs);

            // yy derivative
            eval2.template values<2, false, false>(hessians + n_points, tmp1);
            eval1.template hessians<1, false, false>(tmp1, tmp2);

            // yz derivative
            eval2.template gradients<2, false, false>(hessians + 5 * n_points,
                                                      tmp1);
            eval1.template gradients<1, false, true>(tmp1, tmp2);

            // zz derivative
            eval2.template hessians<2, false, false>(hessians + 2 * n_points,
                                                     tmp1);
            eval1.template values<1, false, true>(tmp1, tmp2);
            eval0.template values<0, false, true>(tmp2, values_dofs);
            break;

          default:
            Assert(false,
                   ExcNotImplemented(
                     "Only 1d, 2d and 3d implemented for Hessian"));
        }
  }



  /**
   * This struct performs the evaluation of function values and gradients for
   * tensor-product finite elements. This is a specialization for elements
   * where the nodal points coincide with the quadrature points like FE_Q
   * shape functions on Gauss-Lobatto elements integrated with Gauss-Lobatto
   * quadrature. The assumption of this class is that the shape 'values'
   * operation is identity, which allows us to write shorter code.
   *
   * In literature, this form of evaluation is often called spectral
   * evaluation, spectral collocation or simply collocation, meaning the same
   * location for shape functions and evaluation space (quadrature points).
   */
  template <int dim, int fe_degree, typename Number>
  struct FEEvaluationImplCollocation
  {
    using Number2 =
      typename FEEvaluationData<dim, Number, false>::shape_info_number_type;
    using Eval = EvaluatorTensorProduct<evaluate_evenodd,
                                        dim,
                                        fe_degree + 1,
                                        fe_degree + 1,
                                        Number,
                                        Number2>;

    static void
    evaluate(const unsigned int                     n_components,
             const EvaluationFlags::EvaluationFlags evaluation_flag,
             const Number                          *values_dofs,
             FEEvaluationData<dim, Number, false>  &fe_eval)
    {
      constexpr std::size_t n_points = Utilities::pow(fe_degree + 1, dim);

      for (unsigned int c = 0; c < n_components; ++c)
        {
          if ((evaluation_flag & EvaluationFlags::values) != 0u)
            for (unsigned int i = 0; i < n_points; ++i)
              fe_eval.begin_values()[n_points * c + i] =
                values_dofs[n_points * c + i];

          if ((evaluation_flag & EvaluationFlags::gradients) != 0u)
            evaluate_gradients_collocation<fe_degree + 1, dim>(
              fe_eval.get_shape_info().data.front(),
              values_dofs + c * n_points,
              fe_eval.begin_gradients() + c * dim * n_points);
        }
    }

    static void
    integrate(const unsigned int                     n_components,
              const EvaluationFlags::EvaluationFlags integration_flag,
              Number                                *values_dofs,
              FEEvaluationData<dim, Number, false>  &fe_eval,
              const bool                             add_into_values_array)
    {
      constexpr std::size_t n_points = Utilities::pow(fe_degree + 1, dim);

      for (unsigned int c = 0; c < n_components; ++c)
        {
          if ((integration_flag & EvaluationFlags::values) != 0u)
            {
              if (add_into_values_array)
                for (unsigned int i = 0; i < n_points; ++i)
                  values_dofs[n_points * c + i] +=
                    fe_eval.begin_values()[n_points * c + i];
              else
                for (unsigned int i = 0; i < n_points; ++i)
                  values_dofs[n_points * c + i] =
                    fe_eval.begin_values()[n_points * c + i];
            }

          if ((integration_flag & EvaluationFlags::gradients) != 0u)
            integrate_gradients_collocation<fe_degree + 1, dim>(
              fe_eval.get_shape_info().data.front(),
              values_dofs + c * n_points,
              fe_eval.begin_gradients() + c * dim * n_points,
              add_into_values_array ||
                ((integration_flag & EvaluationFlags::values) != 0u));
        }
    }
  };



  /**
   * This struct performs the evaluation of function values and gradients for
   * tensor-product finite elements. This is a specialization for symmetric
   * basis functions about the mid point 0.5 of the unit interval with the
   * same number of quadrature points as degrees of freedom. In that case, we
   * can first transform the basis to one that has the nodal points in the
   * quadrature points (i.e., the collocation space) and then perform the
   * evaluation of the first and second derivatives in this transformed space,
   * using the identity operation for the shape values.
   */
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct FEEvaluationImplTransformToCollocation
  {
    static void
    evaluate(const unsigned int                     n_components,
             const EvaluationFlags::EvaluationFlags evaluation_flag,
             const Number                          *values_dofs,
             FEEvaluationData<dim, Number, false>  &fe_eval)
    {
      const auto &shape_data = fe_eval.get_shape_info().data.front();

      Assert(n_q_points_1d > fe_degree,
             ExcMessage("You lose information when going to a collocation "
                        "space of lower degree, so the evaluation results "
                        "would be wrong. Thus, this class does not permit "
                        "the chosen operation."));
      constexpr std::size_t n_dofs     = Utilities::pow(fe_degree + 1, dim);
      constexpr std::size_t n_q_points = Utilities::pow(n_q_points_1d, dim);

      for (unsigned int c = 0; c < n_components; ++c)
        {
          FEEvaluationImplBasisChange<
            evaluate_evenodd,
            EvaluatorQuantity::value,
            dim,
            (fe_degree >= n_q_points_1d ? n_q_points_1d : fe_degree + 1),
            n_q_points_1d>::do_forward(1,
                                       shape_data.shape_values_eo,
                                       values_dofs + c * n_dofs,
                                       fe_eval.begin_values() + c * n_q_points);

          // apply derivatives in the collocation space
          if (evaluation_flag & EvaluationFlags::gradients)
            evaluate_gradients_collocation<n_q_points_1d, dim>(
              shape_data,
              fe_eval.begin_values() + c * n_q_points,
              fe_eval.begin_gradients() + c * dim * n_q_points);
        }
    }

    static void
    integrate(const unsigned int                     n_components,
              const EvaluationFlags::EvaluationFlags integration_flag,
              Number                                *values_dofs,
              FEEvaluationData<dim, Number, false>  &fe_eval,
              const bool                             add_into_values_array)
    {
      const auto &shape_data = fe_eval.get_shape_info().data.front();

      Assert(n_q_points_1d > fe_degree,
             ExcMessage("You lose information when going to a collocation "
                        "space of lower degree, so the evaluation results "
                        "would be wrong. Thus, this class does not permit "
                        "the chosen operation."));
      constexpr std::size_t n_q_points = Utilities::pow(n_q_points_1d, dim);

      for (unsigned int c = 0; c < n_components; ++c)
        {
          // apply derivatives in collocation space
          if (integration_flag & EvaluationFlags::gradients)
            integrate_gradients_collocation<n_q_points_1d, dim>(
              shape_data,
              fe_eval.begin_values() + c * n_q_points,
              fe_eval.begin_gradients() + c * dim * n_q_points,
              /*add_into_values_array=*/
              integration_flag & EvaluationFlags::values);

          // transform back to the original space
          FEEvaluationImplBasisChange<
            evaluate_evenodd,
            EvaluatorQuantity::value,
            dim,
            (fe_degree >= n_q_points_1d ? n_q_points_1d : fe_degree + 1),
            n_q_points_1d>::do_backward(1,
                                        shape_data.shape_values_eo,
                                        add_into_values_array,
                                        fe_eval.begin_values() + c * n_q_points,
                                        values_dofs +
                                          c *
                                            Utilities::pow(fe_degree + 1, dim));
        }
    }
  };



  /**
   * Specialization for MatrixFreeFunctions::tensor_raviart_thomas, which use
   * specific sum-factorization kernels and with normal/tangential shape_data
   */
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct FEEvaluationImpl<MatrixFreeFunctions::tensor_raviart_thomas,
                          dim,
                          fe_degree,
                          n_q_points_1d,
                          Number>
  {
    using Number2 =
      typename FEEvaluationData<dim, Number, false>::shape_info_number_type;

    template <bool integrate>
    static void
    evaluate_or_integrate(
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      Number                                *values_dofs_actual,
      FEEvaluationData<dim, Number, false>  &fe_eval,
      const bool                             add_into_values_array = false);
  };



  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <bool integrate>
  inline void
  FEEvaluationImpl<MatrixFreeFunctions::tensor_raviart_thomas,
                   dim,
                   fe_degree,
                   n_q_points_1d,
                   Number>::
    evaluate_or_integrate(
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      Number                                *values_dofs,
      FEEvaluationData<dim, Number, false>  &fe_eval,
      const bool                             add)
  {
    Assert(dim == 2 || dim == 3,
           ExcMessage("Only dim = 2,3 implemented for Raviart-Thomas "
                      "evaluation/integration"));

    if (evaluation_flag == EvaluationFlags::nothing)
      return;

    AssertDimension(fe_eval.get_shape_info().data.size(), 2);
    AssertDimension(n_q_points_1d,
                    fe_eval.get_shape_info().data[0].n_q_points_1d);
    AssertDimension(n_q_points_1d,
                    fe_eval.get_shape_info().data[1].n_q_points_1d);
    AssertDimension(fe_degree, fe_eval.get_shape_info().data[0].fe_degree);
    AssertDimension(fe_degree, fe_eval.get_shape_info().data[1].fe_degree + 1);

    const auto        &shape_data = fe_eval.get_shape_info().data;
    const unsigned int dofs_per_component =
      Utilities::pow(fe_degree, dim - 1) * (fe_degree + 1);
    const unsigned int n_points  = Utilities::pow(n_q_points_1d, dim);
    Number            *gradients = fe_eval.begin_gradients();
    Number            *values    = fe_eval.begin_values();

    if (integrate)
      {
        EvaluatorTensorProductAnisotropic<dim, fe_degree, n_q_points_1d, false>
          eval;

        const bool do_values = evaluation_flag & EvaluationFlags::values;
        if ((evaluation_flag & EvaluationFlags::gradients) != 0u)
          integrate_gradients_collocation<n_q_points_1d, dim>(shape_data[0],
                                                              values,
                                                              gradients,
                                                              do_values);
        if constexpr (dim > 2)
          eval.template tangential<2, 0>(shape_data[1], values, values);
        eval.template tangential<1, 0>(shape_data[1], values, values);
        eval.template normal<0>(shape_data[0], values, values_dofs, add);

        values += n_points;
        gradients += n_points * dim;
        values_dofs += dofs_per_component;

        if ((evaluation_flag & EvaluationFlags::gradients) != 0u)
          integrate_gradients_collocation<n_q_points_1d, dim>(shape_data[0],
                                                              values,
                                                              gradients,
                                                              do_values);
        if constexpr (dim > 2)
          eval.template tangential<2, 1>(shape_data[1], values, values);
        eval.template tangential<0, 1>(shape_data[1], values, values);
        eval.template normal<1>(shape_data[0], values, values_dofs, add);

        if constexpr (dim > 2)
          {
            values += n_points;
            gradients += n_points * dim;
            values_dofs += dofs_per_component;

            if ((evaluation_flag & EvaluationFlags::gradients) != 0u)
              integrate_gradients_collocation<n_q_points_1d, dim>(shape_data[0],
                                                                  values,
                                                                  gradients,
                                                                  do_values);
            eval.template tangential<1, 2>(shape_data[1], values, values);
            eval.template tangential<0, 2>(shape_data[1], values, values);
            eval.template normal<2>(shape_data[0], values, values_dofs, add);
          }
      }
    else
      {
        EvaluatorTensorProductAnisotropic<dim, fe_degree, n_q_points_1d, true>
          eval;
        eval.template normal<0>(shape_data[0], values_dofs, values);
        eval.template tangential<1, 0>(shape_data[1], values, values);
        if constexpr (dim > 2)
          eval.template tangential<2, 0>(shape_data[1], values, values);
        if ((evaluation_flag & EvaluationFlags::gradients) != 0u)
          evaluate_gradients_collocation<n_q_points_1d, dim>(shape_data[0],
                                                             values,
                                                             gradients);

        values += n_points;
        gradients += n_points * dim;
        values_dofs += dofs_per_component;

        eval.template normal<1>(shape_data[0], values_dofs, values);
        eval.template tangential<0, 1>(shape_data[1], values, values);
        if constexpr (dim > 2)
          eval.template tangential<2, 1>(shape_data[1], values, values);
        if ((evaluation_flag & EvaluationFlags::gradients) != 0u)
          evaluate_gradients_collocation<n_q_points_1d, dim>(shape_data[0],
                                                             values,
                                                             gradients);

        if constexpr (dim > 2)
          {
            values += n_points;
            gradients += n_points * dim;
            values_dofs += dofs_per_component;

            eval.template normal<2>(shape_data[0], values_dofs, values);
            eval.template tangential<0, 2>(shape_data[1], values, values);
            eval.template tangential<1, 2>(shape_data[1], values, values);
            if ((evaluation_flag & EvaluationFlags::gradients) != 0u)
              evaluate_gradients_collocation<n_q_points_1d, dim>(shape_data[0],
                                                                 values,
                                                                 gradients);
          }
      }
  }



  /**
   * This class chooses an appropriate evaluation/integration strategy based on
   * the template parameters and the shape_info variable which contains runtime
   * parameters for the strategy underlying FEEvaluation::evaluate(), i.e.
   * this calls internal::FEEvaluationImpl::evaluate(),
   * internal::FEEvaluationImplCollocation::evaluate() or
   * internal::FEEvaluationImplTransformToCollocation::evaluate() with
   * appropriate template parameters. In case the template parameters
   * fe_degree and n_q_points_1d contain valid information (i.e. fe_degree>-1
   * and n_q_points_1d>0), we simply pass these values to the respective
   * template specializations.  Otherwise, we perform a runtime matching of
   * the runtime parameters to find the correct specialization. This matching
   * currently supports $0\leq fe\_degree \leq 9$ and $degree+1\leq
   * n\_q\_points\_1d\leq fe\_degree+2$.
   */
  template <int dim, typename Number, bool do_integrate>
  struct FEEvaluationImplSelector
  {
    template <int fe_degree, int n_q_points_1d, typename OtherNumber>
    static bool
    run(const unsigned int                     n_components,
        const EvaluationFlags::EvaluationFlags evaluation_flag,
        OtherNumber                           *values_dofs,
        FEEvaluationData<dim, Number, false>  &fe_eval,
        const bool                             sum_into_values_array_in = false)
    {
      // `OtherNumber` is either `const Number` (evaluate()) or `Number`
      // (integrate())
      static_assert(std::is_same_v<Number, std::remove_const_t<OtherNumber>>,
                    "Type of Number and of OtherNumber do not match.");

      const auto element_type = fe_eval.get_shape_info().element_type;
      using ElementType       = MatrixFreeFunctions::ElementType;

      Assert(fe_eval.get_shape_info().data.size() == 1 ||
               (fe_eval.get_shape_info().data.size() == dim &&
                element_type == ElementType::tensor_general) ||
               element_type == ElementType::tensor_raviart_thomas,
             ExcNotImplemented());

      EvaluationFlags::EvaluationFlags actual_flag = evaluation_flag;
      bool sum_into_values_array                   = sum_into_values_array_in;
      if (evaluation_flag & EvaluationFlags::hessians)
        {
          actual_flag |= EvaluationFlags::values;
          Assert(element_type != MatrixFreeFunctions::tensor_none,
                 ExcNotImplemented());
          if constexpr (do_integrate)
            {
              if (fe_eval.get_shape_info().data[0].fe_degree <
                  fe_eval.get_shape_info().data[0].n_q_points_1d)
                integrate_hessians_collocation<n_q_points_1d>(
                  n_components,
                  fe_eval,
                  evaluation_flag & EvaluationFlags::values);
              else
                {
                  integrate_hessians_slow(n_components,
                                          fe_eval,
                                          values_dofs,
                                          sum_into_values_array);
                  sum_into_values_array = true;
                }
            }
        }

      if (fe_degree >= 0 && fe_degree + 1 == n_q_points_1d &&
          element_type == ElementType::tensor_symmetric_collocation)
        {
          evaluate_or_integrate<
            FEEvaluationImplCollocation<dim, fe_degree, Number>>(
            n_components,
            actual_flag,
            values_dofs,
            fe_eval,
            sum_into_values_array);
        }
      // '<=' on type means tensor_symmetric or tensor_symmetric_hermite, see
      // shape_info.h for more details
      else if (fe_degree >= 0 &&
               use_collocation_evaluation(fe_degree, n_q_points_1d) &&
               element_type <= ElementType::tensor_symmetric)
        {
          evaluate_or_integrate<
            FEEvaluationImplTransformToCollocation<dim,
                                                   fe_degree,
                                                   n_q_points_1d,
                                                   Number>>(
            n_components,
            actual_flag,
            values_dofs,
            fe_eval,
            sum_into_values_array);
        }
      else if (fe_degree >= 0 &&
               element_type <= ElementType::tensor_symmetric_no_collocation)
        {
          evaluate_or_integrate<FEEvaluationImpl<ElementType::tensor_symmetric,
                                                 dim,
                                                 fe_degree,
                                                 n_q_points_1d,
                                                 Number>>(
            n_components,
            actual_flag,
            values_dofs,
            fe_eval,
            sum_into_values_array);
        }
      else if (element_type == ElementType::tensor_none)
        {
          evaluate_or_integrate<
            FEEvaluationImpl<ElementType::tensor_none, dim, -1, 0, Number>>(
            n_components,
            actual_flag,
            values_dofs,
            fe_eval,
            sum_into_values_array);
        }
      else if (element_type == ElementType::tensor_symmetric_plus_dg0)
        {
          evaluate_or_integrate<
            FEEvaluationImpl<ElementType::tensor_symmetric_plus_dg0,
                             dim,
                             fe_degree,
                             n_q_points_1d,
                             Number>>(n_components,
                                      actual_flag,
                                      values_dofs,
                                      fe_eval,
                                      sum_into_values_array);
        }
      else if (element_type == ElementType::truncated_tensor)
        {
          evaluate_or_integrate<FEEvaluationImpl<ElementType::truncated_tensor,
                                                 dim,
                                                 fe_degree,
                                                 n_q_points_1d,
                                                 Number>>(
            n_components,
            actual_flag,
            values_dofs,
            fe_eval,
            sum_into_values_array);
        }
      else if (element_type == ElementType::tensor_raviart_thomas)
        {
          if constexpr (fe_degree > 0 && n_q_points_1d > 0 && dim > 1)
            {
              FEEvaluationImpl<ElementType::tensor_raviart_thomas,
                               dim,
                               fe_degree,
                               n_q_points_1d,
                               Number>::
                template evaluate_or_integrate<do_integrate>(
                  actual_flag,
                  const_cast<Number *>(values_dofs),
                  fe_eval,
                  sum_into_values_array);
            }
          else
            {
              Assert(false,
                     ExcNotImplemented("Raviart-Thomas currently only possible "
                                       "in 2d/3d and with templated degree"));
            }
        }
      else
        {
          evaluate_or_integrate<FEEvaluationImpl<ElementType::tensor_general,
                                                 dim,
                                                 fe_degree,
                                                 n_q_points_1d,
                                                 Number>>(
            n_components,
            actual_flag,
            values_dofs,
            fe_eval,
            sum_into_values_array);
        }

      if ((evaluation_flag & EvaluationFlags::hessians) && !do_integrate)
        {
          Assert(element_type != MatrixFreeFunctions::tensor_none,
                 ExcNotImplemented());
          if (fe_eval.get_shape_info().data[0].fe_degree <
              fe_eval.get_shape_info().data[0].n_q_points_1d)
            evaluate_hessians_collocation<n_q_points_1d>(n_components, fe_eval);
          else
            evaluate_hessians_slow(n_components, values_dofs, fe_eval);
        }

      return false;
    }

  private:
    template <typename T>
    static void
    evaluate_or_integrate(
      const unsigned int                     n_components,
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      const Number                          *values_dofs,
      FEEvaluationData<dim, Number, false>  &fe_eval,
      const bool                             sum_into_values_array,
      std::bool_constant<false>)
    {
      (void)sum_into_values_array;

      T::evaluate(n_components, evaluation_flag, values_dofs, fe_eval);
    }

    template <typename T>
    static void
    evaluate_or_integrate(
      const unsigned int                     n_components,
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      Number                                *values_dofs,
      FEEvaluationData<dim, Number, false>  &fe_eval,
      const bool                             sum_into_values_array,
      std::bool_constant<true>)
    {
      T::integrate(n_components,
                   evaluation_flag,
                   values_dofs,
                   fe_eval,
                   sum_into_values_array);
    }

    template <typename T, typename OtherNumber>
    static void
    evaluate_or_integrate(
      const unsigned int                     n_components,
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      OtherNumber                           *values_dofs,
      FEEvaluationData<dim, Number, false>  &fe_eval,
      const bool                             sum_into_values_array)
    {
      evaluate_or_integrate<T>(n_components,
                               evaluation_flag,
                               values_dofs,
                               fe_eval,
                               sum_into_values_array,
                               std::bool_constant<do_integrate>());
    }
  };



  /**
   * This struct implements the action of the inverse @ref GlossMassMatrix "mass matrix" operation,
   * using an FEEvaluationData argument.
   */
  template <int dim, typename Number>
  struct CellwiseInverseMassMatrixImplBasic
  {
    using Number2 =
      typename FEEvaluationData<dim, Number, false>::shape_info_number_type;

    template <int fe_degree, int = 0>
    static bool
    run(const unsigned int                          n_components,
        const FEEvaluationData<dim, Number, false> &fe_eval,
        const Number                               *in_array,
        Number                                     *out_array)
    {
      const unsigned int given_degree =
        (fe_degree > -1) ? fe_degree :
                           fe_eval.get_shape_info().data.front().fe_degree;

      const unsigned int dofs_per_component =
        Utilities::pow(given_degree + 1, dim);

      Assert(dim >= 1 || dim <= 3, ExcNotImplemented());
      Assert(fe_eval.get_shape_info().element_type <=
               MatrixFreeFunctions::tensor_symmetric_no_collocation,
             ExcNotImplemented());

      EvaluatorTensorProduct<evaluate_evenodd,
                             dim,
                             fe_degree + 1,
                             fe_degree + 1,
                             Number,
                             Number2>
        evaluator({},
                  {},
                  fe_eval.get_shape_info().data.front().inverse_shape_values_eo,
                  given_degree + 1,
                  given_degree + 1);

      for (unsigned int d = 0; d < n_components; ++d)
        {
          const Number *in  = in_array + d * dofs_per_component;
          Number       *out = out_array + d * dofs_per_component;
          // Need to select 'apply' method with hessian slot because values
          // assume symmetries that do not exist in the inverse shapes
          evaluator.template hessians<0, true, false>(in, out);
          if (dim > 1)
            evaluator.template hessians<1, true, false>(out, out);
          if (dim > 2)
            evaluator.template hessians<2, true, false>(out, out);
        }
      for (unsigned int q = 0; q < dofs_per_component; ++q)
        {
          const Number inverse_JxW_q = Number(1.) / fe_eval.JxW(q);
          for (unsigned int d = 0; d < n_components; ++d)
            out_array[q + d * dofs_per_component] *= inverse_JxW_q;
        }
      for (unsigned int d = 0; d < n_components; ++d)
        {
          Number *out = out_array + d * dofs_per_component;
          if (dim > 2)
            evaluator.template hessians<2, false, false>(out, out);
          if (dim > 1)
            evaluator.template hessians<1, false, false>(out, out);
          evaluator.template hessians<0, false, false>(out, out);
        }
      return false;
    }
  };



  /**
   * This struct implements the action of the inverse @ref GlossMassMatrix "mass matrix" operation
   * with user-provided coefficients at quadrature points (in contrast to
   * CellwiseInverseMassMatrixImplBasic, which implicitly uses `1/(|J|xW)' as
   * coefficient).
   */
  template <int dim, typename Number>
  struct CellwiseInverseMassMatrixImplFlexible
  {
    using Number2 =
      typename FEEvaluationData<dim, Number, false>::shape_info_number_type;

    template <int fe_degree, int = 0>
    static bool
    run(const unsigned int                          n_desired_components,
        const FEEvaluationData<dim, Number, false> &fe_eval,
        const ArrayView<const Number>              &inverse_coefficients,
        const bool                                  dyadic_coefficients,
        const Number                               *in_array,
        Number                                     *out_array)
    {
      const unsigned int given_degree =
        (fe_degree > -1) ? fe_degree :
                           fe_eval.get_shape_info().data.front().fe_degree;

      const unsigned int dofs_per_component =
        Utilities::pow(given_degree + 1, dim);

      Assert(inverse_coefficients.size() > 0 &&
               inverse_coefficients.size() % dofs_per_component == 0,
             ExcMessage(
               "Expected diagonal to be a multiple of scalar dof per cells"));

      if (!dyadic_coefficients)
        {
          if (inverse_coefficients.size() != dofs_per_component)
            AssertDimension(n_desired_components * dofs_per_component,
                            inverse_coefficients.size());
        }
      else
        {
          AssertDimension(n_desired_components * n_desired_components *
                            dofs_per_component,
                          inverse_coefficients.size());
        }

      Assert(dim >= 1 || dim <= 3, ExcNotImplemented());
      Assert(fe_eval.get_shape_info().element_type <=
               MatrixFreeFunctions::tensor_symmetric_no_collocation,
             ExcNotImplemented());

      EvaluatorTensorProduct<evaluate_evenodd,
                             dim,
                             fe_degree + 1,
                             fe_degree + 1,
                             Number,
                             Number2>
        evaluator({},
                  {},
                  fe_eval.get_shape_info().data.front().inverse_shape_values_eo,
                  given_degree + 1,
                  given_degree + 1);

      const Number *in  = in_array;
      Number       *out = out_array;

      const Number *inv_coefficient = inverse_coefficients.data();

      const unsigned int shift_coefficient =
        inverse_coefficients.size() > dofs_per_component ? dofs_per_component :
                                                           0;

      const auto n_comp_outer = dyadic_coefficients ? 1 : n_desired_components;
      const auto n_comp_inner = dyadic_coefficients ? n_desired_components : 1;

      for (unsigned int d = 0; d < n_comp_outer; ++d)
        {
          for (unsigned int di = 0; di < n_comp_inner; ++di)
            {
              const Number *in_  = in + di * dofs_per_component;
              Number       *out_ = out + di * dofs_per_component;
              evaluator.template hessians<0, true, false>(in_, out_);
              if (dim > 1)
                evaluator.template hessians<1, true, false>(out_, out_);
              if (dim > 2)
                evaluator.template hessians<2, true, false>(out_, out_);
            }
          if (dyadic_coefficients)
            {
              const auto n_coeff_components =
                n_desired_components * n_desired_components;
              if (n_desired_components == dim)
                {
                  for (unsigned int q = 0; q < dofs_per_component; ++q)
                    vmult<dim>(&inv_coefficient[q * n_coeff_components],
                               &in[q],
                               &out[q],
                               dofs_per_component);
                }
              else
                {
                  for (unsigned int q = 0; q < dofs_per_component; ++q)
                    vmult<-1>(&inv_coefficient[q * n_coeff_components],
                              &in[q],
                              &out[q],
                              dofs_per_component,
                              n_desired_components);
                }
            }
          else
            for (unsigned int q = 0; q < dofs_per_component; ++q)
              out[q] *= inv_coefficient[q];

          for (unsigned int di = 0; di < n_comp_inner; ++di)
            {
              Number *out_ = out + di * dofs_per_component;
              if (dim > 2)
                evaluator.template hessians<2, false, false>(out_, out_);
              if (dim > 1)
                evaluator.template hessians<1, false, false>(out_, out_);
              evaluator.template hessians<0, false, false>(out_, out_);
            }

          in += dofs_per_component;
          out += dofs_per_component;
          inv_coefficient += shift_coefficient;
        }

      return false;
    }

  private:
    template <int n_components>
    static inline void
    vmult(const Number      *inverse_coefficients,
          const Number      *src,
          Number            *dst,
          const unsigned int dofs_per_component,
          const unsigned int n_given_components = 0)
    {
      const unsigned int n_desired_components =
        (n_components > -1) ? n_components : n_given_components;

      std::array<Number, dim + 2> tmp = {};
      Assert(n_desired_components <= dim + 2,
             ExcMessage(
               "Number of components larger than dim+2 not supported."));

      for (unsigned int d = 0; d < n_desired_components; ++d)
        tmp[d] = src[d * dofs_per_component];

      for (unsigned int d1 = 0; d1 < n_desired_components; ++d1)
        {
          const Number *inv_coeff_row =
            &inverse_coefficients[d1 * n_desired_components];
          Number sum = inv_coeff_row[0] * tmp[0];
          for (unsigned int d2 = 1; d2 < n_desired_components; ++d2)
            sum += inv_coeff_row[d2] * tmp[d2];
          dst[d1 * dofs_per_component] = sum;
        }
    }
  };



  /**
   * This struct implements the action of a projection of the values given
   * at the quadrature points to the support points,
   * using an FEEvaluationData argument. For the derivation, see comments in
   * step-67.
   */
  template <int dim, typename Number>
  struct CellwiseInverseMassMatrixImplTransformFromQPoints
  {
    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                          n_desired_components,
        const FEEvaluationData<dim, Number, false> &fe_eval,
        const Number                               *in_array,
        Number                                     *out_array)
    {
      static const bool do_inplace =
        fe_degree > -1 && (fe_degree + 1 == n_q_points_1d);

      Assert(fe_eval.get_shape_info().element_type !=
               MatrixFreeFunctions::ElementType::tensor_none,
             ExcNotImplemented());

      const auto &inverse_shape =
        do_inplace ?
          fe_eval.get_shape_info().data.front().inverse_shape_values_eo :
          fe_eval.get_shape_info().data.front().inverse_shape_values;

      const std::size_t dofs_per_component =
        do_inplace ? Utilities::pow(fe_degree + 1, dim) :
                     fe_eval.get_shape_info().dofs_per_component_on_cell;
      const std::size_t n_q_points = do_inplace ?
                                       Utilities::pow(fe_degree + 1, dim) :
                                       fe_eval.get_shape_info().n_q_points;

      using Number2 =
        typename FEEvaluationData<dim, Number, false>::shape_info_number_type;
      EvaluatorTensorProduct<do_inplace ? evaluate_evenodd : evaluate_general,
                             dim,
                             fe_degree + 1,
                             n_q_points_1d,
                             Number,
                             Number2>
        evaluator({},
                  {},
                  inverse_shape,
                  fe_eval.get_shape_info().data.front().fe_degree + 1,
                  fe_eval.get_shape_info().data.front().n_q_points_1d);

      for (unsigned int d = 0; d < n_desired_components; ++d)
        {
          const Number *in  = in_array + d * n_q_points;
          Number       *out = out_array + d * dofs_per_component;

          auto *temp_1 = do_inplace ? out : fe_eval.get_scratch_data().begin();
          auto *temp_2 = do_inplace ?
                           out :
                           (temp_1 + std::max(n_q_points, dofs_per_component));

          if (dim == 3)
            {
              evaluator.template hessians<2, false, false>(in, temp_1);
              evaluator.template hessians<1, false, false>(temp_1, temp_2);
              evaluator.template hessians<0, false, false>(temp_2, out);
            }
          if (dim == 2)
            {
              evaluator.template hessians<1, false, false>(in, temp_1);
              evaluator.template hessians<0, false, false>(temp_1, out);
            }
          if (dim == 1)
            evaluator.template hessians<0, false, false>(in, out);
        }
      return false;
    }
  };
} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
