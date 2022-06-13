// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2022 by the deal.II authors
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


#ifndef dealii_matrix_free_evaluation_kernels_h
#define dealii_matrix_free_evaluation_kernels_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/ndarray.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/matrix_free/dof_info.h>
#include <deal.II/matrix_free/evaluation_flags.h>
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

  template <bool is_long>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_raviart_thomas, is_long>
  {
    static const EvaluatorVariant variant = evaluate_raviart_thomas;
  };



  /**
   * This struct performs the evaluation of function values, gradients and
   * Hessians for tensor-product finite elements. The operation is used for
   * both the symmetric and non-symmetric case, which use different apply
   * functions 'values', 'gradients' in the individual coordinate
   * directions. The apply functions for values are provided through one of
   * the template classes EvaluatorTensorProduct which in turn are selected
   * from the MatrixFreeFunctions::ElementType template argument.
   *
   * There are two specialized implementation classes
   * FEEvaluationImplCollocation (for Gauss-Lobatto elements where the nodal
   * points and the quadrature points coincide and the 'values' operation is
   * identity) and FEEvaluationImplTransformToCollocation (which can be
   * transformed to a collocation space and can then use the identity in these
   * spaces), which both allow for shorter code.
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

    using Eval = EvaluatorTensorProduct<variant,
                                        dim,
                                        fe_degree + 1,
                                        n_q_points_1d,
                                        Number>;

    static void
    evaluate(const unsigned int                     n_components,
             const EvaluationFlags::EvaluationFlags evaluation_flag,
             const Number *                         values_dofs_actual,
             FEEvaluationData<dim, Number, false> & fe_eval);

    static void
    integrate(const unsigned int                     n_components,
              const EvaluationFlags::EvaluationFlags integration_flag,
              Number *                               values_dofs_actual,
              FEEvaluationData<dim, Number, false> & fe_eval,
              const bool                             add_into_values_array);

    static Eval
    create_evaluator_tensor_product(
      const MatrixFreeFunctions::UnivariateShapeData<Number>
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
             const Number *                         values_dofs_actual,
             FEEvaluationData<dim, Number, false> & fe_eval);

    static void
    integrate(const unsigned int                     n_components,
              const EvaluationFlags::EvaluationFlags integration_flag,
              Number *                               values_dofs_actual,
              FEEvaluationData<dim, Number, false> & fe_eval,
              const bool                             add_into_values_array);
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
    template <bool integrate>
    static void
    evaluate_or_integrate(
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      Number *                               values_dofs_actual,
      FEEvaluationData<dim, Number, false> & fe_eval,
      const bool                             add_into_values_array = false);

  private:
    template <typename EvalType>
    static EvalType
    create_evaluator_tensor_product(
      const MatrixFreeFunctions::UnivariateShapeData<Number> &shape_data)
    {
      return EvalType(shape_data.shape_values,
                      shape_data.shape_gradients,
                      shape_data.shape_hessians);
    }

    template <int normal_dir>
    static void
    evaluate_tensor_product_per_component(
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      Number *                               values_dofs_actual,
      FEEvaluationData<dim, Number, false> & fe_eval,
      const bool                             add_into_values_array,
      std::integral_constant<bool, false>);

    template <int normal_dir>
    static void
    evaluate_tensor_product_per_component(
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      Number *                               values_dofs_actual,
      FEEvaluationData<dim, Number, false> & fe_eval,
      const bool                             add_into_values_array,
      std::integral_constant<bool, true>);
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
    const Number *                         values_dofs_actual,
    FEEvaluationData<dim, Number, false> & fe_eval)
  {
    if (evaluation_flag == EvaluationFlags::nothing)
      return;

    std::array<const MatrixFreeFunctions::UnivariateShapeData<Number> *, 3>
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
    const Number *values_dofs = values_dofs_actual;
    if (type == MatrixFreeFunctions::truncated_tensor)
      {
        const std::size_t n_dofs_per_comp =
          fe_eval.get_shape_info().dofs_per_component_on_cell;
        Number *values_dofs_tmp =
          temp1 + 2 * (std::max(n_dofs_per_comp, n_q_points));
        const int degree =
          fe_degree != -1 ? fe_degree : shape_data.front().fe_degree;
        for (unsigned int c = 0; c < n_components; ++c)
          for (int i = 0, count_p = 0, count_q = 0;
               i < (dim > 2 ? degree + 1 : 1);
               ++i)
            {
              for (int j = 0; j < (dim > 1 ? degree + 1 - i : 1); ++j)
                {
                  for (int k = 0; k < degree + 1 - j - i;
                       ++k, ++count_p, ++count_q)
                    values_dofs_tmp[c * dofs_per_comp + count_q] =
                      values_dofs_actual[c * n_dofs_per_comp + count_p];
                  for (int k = degree + 1 - j - i; k < degree + 1;
                       ++k, ++count_q)
                    values_dofs_tmp[c * dofs_per_comp + count_q] = Number();
                }
              for (int j = degree + 1 - i; j < degree + 1; ++j)
                for (int k = 0; k < degree + 1; ++k, ++count_q)
                  values_dofs_tmp[c * dofs_per_comp + count_q] = Number();
            }
        values_dofs = values_dofs_tmp;
      }

    Number *values_quad    = fe_eval.begin_values();
    Number *gradients_quad = fe_eval.begin_gradients();
    Number *hessians_quad  = fe_eval.begin_hessians();

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
              if (evaluation_flag & EvaluationFlags::hessians)
                eval0.template hessians<0, true, false>(values_dofs,
                                                        hessians_quad);

              // advance the next component in 1D array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += n_q_points;
              hessians_quad += n_q_points;
            }
          break;

        case 2:
          for (unsigned int c = 0; c < n_components; ++c)
            {
              // grad x
              if (evaluation_flag & EvaluationFlags::gradients)
                {
                  eval0.template gradients<0, true, false>(values_dofs, temp1);
                  eval1.template values<1, true, false>(temp1, gradients_quad);
                }
              if (evaluation_flag & EvaluationFlags::hessians)
                {
                  // grad xy
                  if (!(evaluation_flag & EvaluationFlags::gradients))
                    eval0.template gradients<0, true, false>(values_dofs,
                                                             temp1);
                  eval1.template gradients<1, true, false>(temp1,
                                                           hessians_quad +
                                                             2 * n_q_points);

                  // grad xx
                  eval0.template hessians<0, true, false>(values_dofs, temp1);
                  eval1.template values<1, true, false>(temp1, hessians_quad);
                }

              // grad y
              eval0.template values<0, true, false>(values_dofs, temp1);
              if (evaluation_flag & EvaluationFlags::gradients)
                eval1.template gradients<1, true, false>(temp1,
                                                         gradients_quad +
                                                           n_q_points);

              // grad yy
              if (evaluation_flag & EvaluationFlags::hessians)
                eval1.template hessians<1, true, false>(temp1,
                                                        hessians_quad +
                                                          n_q_points);

              // val: can use values applied in x
              if (evaluation_flag & EvaluationFlags::values)
                eval1.template values<1, true, false>(temp1, values_quad);

              // advance to the next component in 1D array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += 2 * n_q_points;
              hessians_quad += 3 * n_q_points;
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
                  eval2.template values<2, true, false>(temp2, gradients_quad);
                }

              if (evaluation_flag & EvaluationFlags::hessians)
                {
                  // grad xz
                  if (!(evaluation_flag & EvaluationFlags::gradients))
                    {
                      eval0.template gradients<0, true, false>(values_dofs,
                                                               temp1);
                      eval1.template values<1, true, false>(temp1, temp2);
                    }
                  eval2.template gradients<2, true, false>(temp2,
                                                           hessians_quad +
                                                             4 * n_q_points);

                  // grad xy
                  eval1.template gradients<1, true, false>(temp1, temp2);
                  eval2.template values<2, true, false>(temp2,
                                                        hessians_quad +
                                                          3 * n_q_points);

                  // grad xx
                  eval0.template hessians<0, true, false>(values_dofs, temp1);
                  eval1.template values<1, true, false>(temp1, temp2);
                  eval2.template values<2, true, false>(temp2, hessians_quad);
                }

              // grad y
              eval0.template values<0, true, false>(values_dofs, temp1);
              if (evaluation_flag & EvaluationFlags::gradients)
                {
                  eval1.template gradients<1, true, false>(temp1, temp2);
                  eval2.template values<2, true, false>(temp2,
                                                        gradients_quad +
                                                          n_q_points);
                }

              if (evaluation_flag & EvaluationFlags::hessians)
                {
                  // grad yz
                  if (!(evaluation_flag & EvaluationFlags::gradients))
                    eval1.template gradients<1, true, false>(temp1, temp2);
                  eval2.template gradients<2, true, false>(temp2,
                                                           hessians_quad +
                                                             5 * n_q_points);

                  // grad yy
                  eval1.template hessians<1, true, false>(temp1, temp2);
                  eval2.template values<2, true, false>(temp2,
                                                        hessians_quad +
                                                          n_q_points);
                }

              // grad z: can use the values applied in x direction stored in
              // temp1
              eval1.template values<1, true, false>(temp1, temp2);
              if (evaluation_flag & EvaluationFlags::gradients)
                eval2.template gradients<2, true, false>(temp2,
                                                         gradients_quad +
                                                           2 * n_q_points);

              // grad zz: can use the values applied in x and y direction stored
              // in temp2
              if (evaluation_flag & EvaluationFlags::hessians)
                eval2.template hessians<2, true, false>(temp2,
                                                        hessians_quad +
                                                          2 * n_q_points);

              // val: can use the values applied in x & y direction stored in
              // temp2
              if (evaluation_flag & EvaluationFlags::values)
                eval2.template values<2, true, false>(temp2, values_quad);

              // advance to the next component in 1D array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += 3 * n_q_points;
              hessians_quad += 6 * n_q_points;
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
    Number *                               values_dofs_actual,
    FEEvaluationData<dim, Number, false> & fe_eval,
    const bool                             add_into_values_array)
  {
    std::array<const MatrixFreeFunctions::UnivariateShapeData<Number> *, 3>
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
    Number *hessians_quad  = fe_eval.begin_hessians();

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
              if ((integration_flag & EvaluationFlags::hessians) != 0u)
                {
                  if ((integration_flag & EvaluationFlags::values) != 0u ||
                      (integration_flag & EvaluationFlags::gradients) != 0u ||
                      add_into_values_array == true)
                    eval0.template hessians<0, false, true>(hessians_quad,
                                                            values_dofs);
                  else
                    eval0.template hessians<0, false, false>(hessians_quad,
                                                             values_dofs);
                }

              // advance to the next component in 1D array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += n_q_points;
              hessians_quad += n_q_points;
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
                  eval1.template gradients<1, false, false>(gradients_quad +
                                                              n_q_points,
                                                            temp1);
                  if (integration_flag & EvaluationFlags::values)
                    eval1.template values<1, false, true>(values_quad, temp1);
                  if (add_into_values_array == false)
                    eval0.template values<0, false, false>(temp1, values_dofs);
                  else
                    eval0.template values<0, false, true>(temp1, values_dofs);
                  eval1.template values<1, false, false>(gradients_quad, temp1);
                  eval0.template gradients<0, false, true>(temp1, values_dofs);
                }
              if ((integration_flag & EvaluationFlags::hessians) != 0u)
                {
                  // grad xx
                  eval1.template values<1, false, false>(hessians_quad, temp1);

                  if ((integration_flag & EvaluationFlags::values) != 0u ||
                      (integration_flag & EvaluationFlags::gradients) != 0u ||
                      add_into_values_array == true)
                    eval0.template hessians<0, false, true>(temp1, values_dofs);
                  else
                    eval0.template hessians<0, false, false>(temp1,
                                                             values_dofs);

                  // grad yy
                  eval1.template hessians<1, false, false>(hessians_quad +
                                                             n_q_points,
                                                           temp1);
                  eval0.template values<0, false, true>(temp1, values_dofs);

                  // grad xy
                  eval1.template gradients<1, false, false>(hessians_quad +
                                                              2 * n_q_points,
                                                            temp1);
                  eval0.template gradients<0, false, true>(temp1, values_dofs);
                }

              // advance to the next component in 1D array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += 2 * n_q_points;
              hessians_quad += 3 * n_q_points;
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
                  eval2.template gradients<2, false, false>(gradients_quad +
                                                              2 * n_q_points,
                                                            temp1);
                  if (integration_flag & EvaluationFlags::values)
                    eval2.template values<2, false, true>(values_quad, temp1);
                  eval1.template values<1, false, false>(temp1, temp2);
                  eval2.template values<2, false, false>(gradients_quad +
                                                           n_q_points,
                                                         temp1);
                  eval1.template gradients<1, false, true>(temp1, temp2);
                  if (add_into_values_array == false)
                    eval0.template values<0, false, false>(temp2, values_dofs);
                  else
                    eval0.template values<0, false, true>(temp2, values_dofs);
                  eval2.template values<2, false, false>(gradients_quad, temp1);
                  eval1.template values<1, false, false>(temp1, temp2);
                  eval0.template gradients<0, false, true>(temp2, values_dofs);
                }
              if ((integration_flag & EvaluationFlags::hessians) != 0u)
                {
                  // grad xx
                  eval2.template values<2, false, false>(hessians_quad, temp1);
                  eval1.template values<1, false, false>(temp1, temp2);

                  if ((integration_flag & EvaluationFlags::values) != 0u ||
                      (integration_flag & EvaluationFlags::gradients) != 0u ||
                      add_into_values_array == true)
                    eval0.template hessians<0, false, true>(temp2, values_dofs);
                  else
                    eval0.template hessians<0, false, false>(temp2,
                                                             values_dofs);

                  // grad yy
                  eval2.template values<2, false, false>(hessians_quad +
                                                           n_q_points,
                                                         temp1);
                  eval1.template hessians<1, false, false>(temp1, temp2);
                  eval0.template values<0, false, true>(temp2, values_dofs);

                  // grad zz
                  eval2.template hessians<2, false, false>(hessians_quad +
                                                             2 * n_q_points,
                                                           temp1);
                  eval1.template values<1, false, false>(temp1, temp2);
                  eval0.template values<0, false, true>(temp2, values_dofs);

                  // grad xy
                  eval2.template values<2, false, false>(hessians_quad +
                                                           3 * n_q_points,
                                                         temp1);
                  eval1.template gradients<1, false, false>(temp1, temp2);
                  eval0.template gradients<0, false, true>(temp2, values_dofs);

                  // grad xz
                  eval2.template gradients<2, false, false>(hessians_quad +
                                                              4 * n_q_points,
                                                            temp1);
                  eval1.template values<1, false, false>(temp1, temp2);
                  eval0.template gradients<0, false, true>(temp2, values_dofs);

                  // grad yz
                  eval2.template gradients<2, false, false>(hessians_quad +
                                                              5 * n_q_points,
                                                            temp1);
                  eval1.template gradients<1, false, false>(temp1, temp2);
                  eval0.template values<0, false, true>(temp2, values_dofs);
                }

              // advance to the next component in 1D array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += 3 * n_q_points;
              hessians_quad += 6 * n_q_points;
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
      {
        const std::size_t n_dofs_per_comp =
          fe_eval.get_shape_info().dofs_per_component_on_cell;
        values_dofs -= dofs_per_comp * n_components;
        const int degree =
          fe_degree != -1 ? fe_degree : shape_data.front().fe_degree;
        for (unsigned int c = 0; c < n_components; ++c)
          for (int i = 0, count_p = 0, count_q = 0;
               i < (dim > 2 ? degree + 1 : 1);
               ++i)
            {
              for (int j = 0; j < (dim > 1 ? degree + 1 - i : 1); ++j)
                {
                  for (int k = 0; k < degree + 1 - j - i;
                       ++k, ++count_p, ++count_q)
                    values_dofs_actual[c * n_dofs_per_comp + count_p] =
                      values_dofs[c * dofs_per_comp + count_q];
                  count_q += j + i;
                }
              count_q += i * (degree + 1);
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
    Number>::evaluate(const unsigned int                     n_components,
                      const EvaluationFlags::EvaluationFlags evaluation_flag,
                      const Number *                         values_dofs_actual,
                      FEEvaluationData<dim, Number, false> & fe_eval)
  {
    const std::size_t n_dofs =
      fe_eval.get_shape_info().dofs_per_component_on_cell;
    const std::size_t n_q_points = fe_eval.get_shape_info().n_q_points;

    const auto &shape_data = fe_eval.get_shape_info().data;

    using Eval =
      EvaluatorTensorProduct<evaluate_general, 1, 0, 0, Number, Number>;

    if (evaluation_flag & EvaluationFlags::values)
      {
        const auto shape_values    = shape_data.front().shape_values.data();
        auto       values_quad_ptr = fe_eval.begin_values();
        auto       values_dofs_actual_ptr = values_dofs_actual;

        Eval eval(shape_values, nullptr, nullptr, n_dofs, n_q_points);
        for (unsigned int c = 0; c < n_components; ++c)
          {
            eval.template values<0, true, false>(values_dofs_actual_ptr,
                                                 values_quad_ptr);

            values_quad_ptr += n_q_points;
            values_dofs_actual_ptr += n_dofs;
          }
      }

    if (evaluation_flag & EvaluationFlags::gradients)
      {
        const auto shape_gradients = shape_data.front().shape_gradients.data();
        auto       gradients_quad_ptr     = fe_eval.begin_gradients();
        auto       values_dofs_actual_ptr = values_dofs_actual;

        for (unsigned int c = 0; c < n_components; ++c)
          {
            for (unsigned int d = 0; d < dim; ++d)
              {
                Eval eval(nullptr,
                          shape_gradients + n_q_points * n_dofs * d,
                          nullptr,
                          n_dofs,
                          n_q_points);

                eval.template gradients<0, true, false>(values_dofs_actual_ptr,
                                                        gradients_quad_ptr);

                gradients_quad_ptr += n_q_points;
              }
            values_dofs_actual_ptr += n_dofs;
          }
      }

    if (evaluation_flag & EvaluationFlags::hessians)
      Assert(false, ExcNotImplemented());
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
                       Number *                              values_dofs_actual,
                       FEEvaluationData<dim, Number, false> &fe_eval,
                       const bool add_into_values_array)
  {
    // TODO: implement hessians
    AssertThrow(!(integration_flag & EvaluationFlags::hessians),
                ExcNotImplemented());

    const std::size_t n_dofs =
      fe_eval.get_shape_info().dofs_per_component_on_cell;
    const std::size_t n_q_points = fe_eval.get_shape_info().n_q_points;

    const auto &shape_data = fe_eval.get_shape_info().data;

    using Eval =
      EvaluatorTensorProduct<evaluate_general, 1, 0, 0, Number, Number>;

    if (integration_flag & EvaluationFlags::values)
      {
        const auto shape_values    = shape_data.front().shape_values.data();
        auto       values_quad_ptr = fe_eval.begin_values();
        auto       values_dofs_actual_ptr = values_dofs_actual;

        Eval eval(shape_values, nullptr, nullptr, n_dofs, n_q_points);
        for (unsigned int c = 0; c < n_components; ++c)
          {
            if (add_into_values_array == false)
              eval.template values<0, false, false>(values_quad_ptr,
                                                    values_dofs_actual_ptr);
            else
              eval.template values<0, false, true>(values_quad_ptr,
                                                   values_dofs_actual_ptr);

            values_quad_ptr += n_q_points;
            values_dofs_actual_ptr += n_dofs;
          }
      }

    if (integration_flag & EvaluationFlags::gradients)
      {
        const auto shape_gradients = shape_data.front().shape_gradients.data();
        auto       gradients_quad_ptr     = fe_eval.begin_gradients();
        auto       values_dofs_actual_ptr = values_dofs_actual;

        for (unsigned int c = 0; c < n_components; ++c)
          {
            for (unsigned int d = 0; d < dim; ++d)
              {
                Eval eval(nullptr,
                          shape_gradients + n_q_points * n_dofs * d,
                          nullptr,
                          n_dofs,
                          n_q_points);

                if ((add_into_values_array == false &&
                     !(integration_flag & EvaluationFlags::values)) &&
                    d == 0)
                  eval.template gradients<0, false, false>(
                    gradients_quad_ptr, values_dofs_actual_ptr);
                else
                  eval.template gradients<0, false, true>(
                    gradients_quad_ptr, values_dofs_actual_ptr);

                gradients_quad_ptr += n_q_points;
              }
            values_dofs_actual_ptr += n_dofs;
          }
      }
  }


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
      Number *                               values_dofs_actual,
      FEEvaluationData<dim, Number, false> & fe_eval,
      const bool                             add_into_values_array)
  {
    if (evaluation_flag == EvaluationFlags::nothing)
      return;

    AssertDimension(fe_eval.get_shape_info().data.size(), 2);
    // First component:
    evaluate_tensor_product_per_component<0>(
      evaluation_flag,
      values_dofs_actual,
      fe_eval,
      add_into_values_array,
      std::integral_constant<bool, integrate>());
    // Second component :
    evaluate_tensor_product_per_component<1>(
      evaluation_flag,
      values_dofs_actual,
      fe_eval,
      add_into_values_array,
      std::integral_constant<bool, integrate>());
    if (dim == 3)
      {
        // Third component
        evaluate_tensor_product_per_component<2>(
          evaluation_flag,
          values_dofs_actual,
          fe_eval,
          add_into_values_array,
          std::integral_constant<bool, integrate>());
      }
  }

  // Helper function that applies the 1d evaluation kernels.
  // std::integral_constant<bool, false> is the interpolation path, and
  // std::integral_constant<bool, true> below is the integration path.
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int normal_dir>
  inline void
  FEEvaluationImpl<MatrixFreeFunctions::tensor_raviart_thomas,
                   dim,
                   fe_degree,
                   n_q_points_1d,
                   Number>::
    evaluate_tensor_product_per_component(
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      Number *                               values_dofs_actual,
      FEEvaluationData<dim, Number, false> & fe_eval,
      const bool                             add_into_values_array,
      std::integral_constant<bool, false>)
  {
    (void)add_into_values_array;

    using EvalNormal =
      EvaluatorTensorProductAnisotropic<evaluate_raviart_thomas,
                                        dim,
                                        (fe_degree == -1) ? 1 : fe_degree + 1,
                                        n_q_points_1d,
                                        Number,
                                        normal_dir>;

    using EvalTangent =
      EvaluatorTensorProductAnisotropic<evaluate_raviart_thomas,
                                        dim,
                                        (fe_degree == -1) ? 1 : fe_degree,
                                        n_q_points_1d,
                                        Number,
                                        normal_dir>;
    using Eval0 =
      typename std::conditional<normal_dir == 0, EvalNormal, EvalTangent>::type;
    using Eval1 =
      typename std::conditional<normal_dir == 1, EvalNormal, EvalTangent>::type;
    using Eval2 =
      typename std::conditional<normal_dir == 2, EvalNormal, EvalTangent>::type;

    const MatrixFreeFunctions::ShapeInfo<Number> &shape_info =
      fe_eval.get_shape_info();
    Eval0 eval0 = create_evaluator_tensor_product<Eval0>(
      ((normal_dir == 0) ? shape_info.data[0] : shape_info.data[1]));
    Eval1 eval1 = create_evaluator_tensor_product<Eval1>(
      ((normal_dir == 1) ? shape_info.data[0] : shape_info.data[1]));
    Eval2 eval2 = create_evaluator_tensor_product<Eval2>(
      ((normal_dir == 2) ? shape_info.data[0] : shape_info.data[1]));

    Number *temp1 = fe_eval.get_scratch_data().begin();
    Number *temp2;

    temp2 =
      temp1 +
      std::max(Utilities::fixed_power<dim>(shape_info.data[0].fe_degree + 1),
               Utilities::fixed_power<dim>(shape_info.data[0].n_q_points_1d));

    const std::size_t n_q_points    = shape_info.n_q_points;
    const std::size_t dofs_per_comp = shape_info.dofs_per_component_on_cell;

    // Initial shift depending on component (normal_dir)
    Number *values_dofs = values_dofs_actual + dofs_per_comp * normal_dir;
    Number *values_quad = fe_eval.begin_values() + n_q_points * normal_dir;
    Number *gradients_quad =
      fe_eval.begin_gradients() + dim * n_q_points * normal_dir;
    Number *hessians_quad =
      (dim == 2) ? fe_eval.begin_hessians() + 3 * n_q_points * normal_dir :
                   fe_eval.begin_hessians() + 6 * n_q_points * normal_dir;

    switch (dim)
      {
        case 2:
          if (evaluation_flag & EvaluationFlags::gradients)
            {
              eval0.template gradients<0, true, false>(values_dofs, temp1);
              eval1.template values<1, true, false>(temp1, gradients_quad);
            }
          if (evaluation_flag & EvaluationFlags::hessians)
            {
              // The evaluation/integration here *should* work, however
              // the piola transform is not implemented.
              AssertThrow(false, ExcNotImplemented());
              // grad xy
              if (!(evaluation_flag & EvaluationFlags::gradients))
                eval0.template gradients<0, true, false>(values_dofs, temp1);
              eval1.template gradients<1, true, false>(temp1,
                                                       hessians_quad +
                                                         2 * n_q_points);

              // grad xx
              eval0.template hessians<0, true, false>(values_dofs, temp1);
              eval1.template values<1, true, false>(temp1, hessians_quad);
            }

          // grad y
          eval0.template values<0, true, false>(values_dofs, temp1);
          if (evaluation_flag & EvaluationFlags::gradients)
            eval1.template gradients<1, true, false>(temp1,
                                                     gradients_quad +
                                                       n_q_points);

          // grad yy
          if (evaluation_flag & EvaluationFlags::hessians)
            eval1.template hessians<1, true, false>(temp1,
                                                    hessians_quad + n_q_points);

          // val: can use values applied in x
          if (evaluation_flag & EvaluationFlags::values)
            eval1.template values<1, true, false>(temp1, values_quad);
          break;
        case 3:
          if (evaluation_flag & EvaluationFlags::gradients)
            {
              // grad x
              eval0.template gradients<0, true, false>(values_dofs, temp1);
              eval1.template values<1, true, false>(temp1, temp2);
              eval2.template values<2, true, false>(temp2, gradients_quad);
            }

          if (evaluation_flag & EvaluationFlags::hessians)
            {
              // The evaluation/integration here *should* work, however
              // the piola transform is not implemented.
              AssertThrow(false, ExcNotImplemented());
              // grad xz
              if (!(evaluation_flag & EvaluationFlags::gradients))
                {
                  eval0.template gradients<0, true, false>(values_dofs, temp1);
                  eval1.template values<1, true, false>(temp1, temp2);
                }
              eval2.template gradients<2, true, false>(temp2,
                                                       hessians_quad +
                                                         4 * n_q_points);

              // grad xy
              eval1.template gradients<1, true, false>(temp1, temp2);
              eval2.template values<2, true, false>(temp2,
                                                    hessians_quad +
                                                      3 * n_q_points);

              // grad xx
              eval0.template hessians<0, true, false>(values_dofs, temp1);
              eval1.template values<1, true, false>(temp1, temp2);
              eval2.template values<2, true, false>(temp2, hessians_quad);
            }

          // grad y
          eval0.template values<0, true, false>(values_dofs, temp1);
          if (evaluation_flag & EvaluationFlags::gradients)
            {
              eval1.template gradients<1, true, false>(temp1, temp2);
              eval2.template values<2, true, false>(temp2,
                                                    gradients_quad +
                                                      n_q_points);
            }

          if (evaluation_flag & EvaluationFlags::hessians)
            {
              // grad yz
              if (!(evaluation_flag & EvaluationFlags::gradients))
                eval1.template gradients<1, true, false>(temp1, temp2);
              eval2.template gradients<2, true, false>(temp2,
                                                       hessians_quad +
                                                         5 * n_q_points);

              // grad yy
              eval1.template hessians<1, true, false>(temp1, temp2);
              eval2.template values<2, true, false>(temp2,
                                                    hessians_quad + n_q_points);
            }

          // grad z: can use the values applied in x direction stored in
          // temp1
          eval1.template values<1, true, false>(temp1, temp2);
          if (evaluation_flag & EvaluationFlags::gradients)
            eval2.template gradients<2, true, false>(temp2,
                                                     gradients_quad +
                                                       2 * n_q_points);

          // grad zz: can use the values applied in x and y direction stored
          // in temp2
          if (evaluation_flag & EvaluationFlags::hessians)
            eval2.template hessians<2, true, false>(temp2,
                                                    hessians_quad +
                                                      2 * n_q_points);

          // val: can use the values applied in x & y direction stored in
          // temp2
          if (evaluation_flag & EvaluationFlags::values)
            eval2.template values<2, true, false>(temp2, values_quad);
          break;
        default:
          AssertThrow(false, ExcNotImplemented());
      }
  }

  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int normal_dir>
  inline void
  FEEvaluationImpl<MatrixFreeFunctions::tensor_raviart_thomas,
                   dim,
                   fe_degree,
                   n_q_points_1d,
                   Number>::
    evaluate_tensor_product_per_component(
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      Number *                               values_dofs_actual,
      FEEvaluationData<dim, Number, false> & fe_eval,
      const bool                             add_into_values_array,
      std::integral_constant<bool, true>)
  {
    using EvalNormal =
      EvaluatorTensorProductAnisotropic<evaluate_raviart_thomas,
                                        dim,
                                        (fe_degree == -1) ? 1 : fe_degree + 1,
                                        n_q_points_1d,
                                        Number,
                                        normal_dir>;

    using EvalTangent =
      EvaluatorTensorProductAnisotropic<evaluate_raviart_thomas,
                                        dim,
                                        (fe_degree == -1) ? 1 : fe_degree,
                                        n_q_points_1d,
                                        Number,
                                        normal_dir>;
    using Eval0 =
      typename std::conditional<normal_dir == 0, EvalNormal, EvalTangent>::type;
    using Eval1 =
      typename std::conditional<normal_dir == 1, EvalNormal, EvalTangent>::type;
    using Eval2 =
      typename std::conditional<normal_dir == 2, EvalNormal, EvalTangent>::type;

    const MatrixFreeFunctions::ShapeInfo<Number> &shape_info =
      fe_eval.get_shape_info();
    Eval0 eval0 = create_evaluator_tensor_product<Eval0>(
      ((normal_dir == 0) ? shape_info.data[0] : shape_info.data[1]));
    Eval1 eval1 = create_evaluator_tensor_product<Eval1>(
      ((normal_dir == 1) ? shape_info.data[0] : shape_info.data[1]));
    Eval2 eval2 = create_evaluator_tensor_product<Eval2>(
      ((normal_dir == 2) ? shape_info.data[0] : shape_info.data[1]));

    Number *temp1 = fe_eval.get_scratch_data().begin();
    Number *temp2;

    temp2 =
      temp1 +
      std::max(Utilities::fixed_power<dim>(shape_info.data[0].fe_degree + 1),
               Utilities::fixed_power<dim>(shape_info.data[0].n_q_points_1d));

    const std::size_t n_q_points    = shape_info.n_q_points;
    const std::size_t dofs_per_comp = shape_info.dofs_per_component_on_cell;

    // Initial shift depending on component (normal_dir)
    Number *values_dofs = values_dofs_actual + dofs_per_comp * normal_dir;
    Number *values_quad = fe_eval.begin_values() + n_q_points * normal_dir;
    Number *gradients_quad =
      fe_eval.begin_gradients() + dim * n_q_points * normal_dir;
    Number *hessians_quad =
      (dim == 2) ? fe_eval.begin_hessians() + 3 * n_q_points * normal_dir :
                   fe_eval.begin_hessians() + 6 * n_q_points * normal_dir;

    // Integrate path
    switch (dim)
      {
        case 2:
          if ((evaluation_flag & EvaluationFlags::values) &&
              !(evaluation_flag & EvaluationFlags::gradients))
            {
              eval1.template values<1, false, false>(values_quad, temp1);
              if (add_into_values_array == false)
                eval0.template values<0, false, false>(temp1, values_dofs);
              else
                eval0.template values<0, false, true>(temp1, values_dofs);
            }
          if (evaluation_flag & EvaluationFlags::gradients)
            {
              eval1.template gradients<1, false, false>(gradients_quad +
                                                          n_q_points,
                                                        temp1);
              if ((evaluation_flag & EvaluationFlags::values))
                eval1.template values<1, false, true>(values_quad, temp1);
              if (add_into_values_array == false)
                eval0.template values<0, false, false>(temp1, values_dofs);
              else
                eval0.template values<0, false, true>(temp1, values_dofs);
              eval1.template values<1, false, false>(gradients_quad, temp1);
              eval0.template gradients<0, false, true>(temp1, values_dofs);
            }
          if (evaluation_flag & EvaluationFlags::hessians)
            {
              // grad xx
              eval1.template values<1, false, false>(hessians_quad, temp1);

              if ((evaluation_flag & EvaluationFlags::values) ||
                  (evaluation_flag & EvaluationFlags::gradients) ||
                  add_into_values_array == true)
                eval0.template hessians<0, false, true>(temp1, values_dofs);
              else
                eval0.template hessians<0, false, false>(temp1, values_dofs);

              // grad yy
              eval1.template hessians<1, false, false>(hessians_quad +
                                                         n_q_points,
                                                       temp1);
              eval0.template values<0, false, true>(temp1, values_dofs);

              // grad xy
              eval1.template gradients<1, false, false>(hessians_quad +
                                                          2 * n_q_points,
                                                        temp1);
              eval0.template gradients<0, false, true>(temp1, values_dofs);
            }
          break;

        case 3:
          if ((evaluation_flag & EvaluationFlags::values) &&
              !(evaluation_flag & EvaluationFlags::gradients))
            {
              eval2.template values<2, false, false>(values_quad, temp1);
              eval1.template values<1, false, false>(temp1, temp2);
              if (add_into_values_array == false)
                eval0.template values<0, false, false>(temp2, values_dofs);
              else
                eval0.template values<0, false, true>(temp2, values_dofs);
            }
          if (evaluation_flag & EvaluationFlags::gradients)
            {
              eval2.template gradients<2, false, false>(gradients_quad +
                                                          2 * n_q_points,
                                                        temp1);
              if ((evaluation_flag & EvaluationFlags::values))
                eval2.template values<2, false, true>(values_quad, temp1);
              eval1.template values<1, false, false>(temp1, temp2);
              eval2.template values<2, false, false>(gradients_quad +
                                                       n_q_points,
                                                     temp1);
              eval1.template gradients<1, false, true>(temp1, temp2);
              if (add_into_values_array == false)
                eval0.template values<0, false, false>(temp2, values_dofs);
              else
                eval0.template values<0, false, true>(temp2, values_dofs);
              eval2.template values<2, false, false>(gradients_quad, temp1);
              eval1.template values<1, false, false>(temp1, temp2);
              eval0.template gradients<0, false, true>(temp2, values_dofs);
            }
          if (evaluation_flag & EvaluationFlags::hessians)
            {
              // grad xx
              eval2.template values<2, false, false>(hessians_quad, temp1);
              eval1.template values<1, false, false>(temp1, temp2);

              if ((evaluation_flag & EvaluationFlags::values) ||
                  (evaluation_flag & EvaluationFlags::gradients) ||
                  add_into_values_array == true)
                eval0.template hessians<0, false, true>(temp2, values_dofs);
              else
                eval0.template hessians<0, false, false>(temp2, values_dofs);

              // grad yy
              eval2.template values<2, false, false>(hessians_quad + n_q_points,
                                                     temp1);
              eval1.template hessians<1, false, false>(temp1, temp2);
              eval0.template values<0, false, true>(temp2, values_dofs);

              // grad zz
              eval2.template hessians<2, false, false>(hessians_quad +
                                                         2 * n_q_points,
                                                       temp1);
              eval1.template values<1, false, false>(temp1, temp2);
              eval0.template values<0, false, true>(temp2, values_dofs);

              // grad xy
              eval2.template values<2, false, false>(hessians_quad +
                                                       3 * n_q_points,
                                                     temp1);
              eval1.template gradients<1, false, false>(temp1, temp2);
              eval0.template gradients<0, false, true>(temp2, values_dofs);

              // grad xz
              eval2.template gradients<2, false, false>(hessians_quad +
                                                          4 * n_q_points,
                                                        temp1);
              eval1.template values<1, false, false>(temp1, temp2);
              eval0.template gradients<0, false, true>(temp2, values_dofs);

              // grad yz
              eval2.template gradients<2, false, false>(hessians_quad +
                                                          5 * n_q_points,
                                                        temp1);
              eval1.template gradients<1, false, false>(temp1, temp2);
              eval0.template values<0, false, true>(temp2, values_dofs);
            }

          break;
        default:
          AssertThrow(false, ExcNotImplemented());
      }
  }

  /**
   * This struct implements the change between two different bases. This is an
   * ingredient in the FEEvaluationImplTransformToCollocation class where we
   * first transform to the appropriate basis where we can compute the
   * derivative through collocation techniques.
   *
   * This class allows for dimension-independent application of the operation,
   * implemented by template recursion. It has been tested up to 6D.
   */
  template <EvaluatorVariant  variant,
            EvaluatorQuantity quantity,
            int               dim,
            int               basis_size_1,
            int               basis_size_2,
            typename Number,
            typename Number2>
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
     * @param basis_size_1_variable In case the template argument basis_size_1
     * is zero, the size of the first basis can alternatively be passed in as a
     * run time argument. The template argument takes precedence in case it is
     * nonzero for efficiency reasons.
     * @param basis_size_2_variable In case the template argument basis_size_1
     * is zero, the size of the second basis can alternatively be passed in as a
     * run time argument.
     */
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
    static void
    do_forward(
      const unsigned int            n_components,
      const AlignedVector<Number2> &transformation_matrix,
      const Number *                values_in,
      Number *                      values_out,
      const unsigned int basis_size_1_variable = numbers::invalid_unsigned_int,
      const unsigned int basis_size_2_variable = numbers::invalid_unsigned_int)
    {
      Assert(
        basis_size_1 != 0 || basis_size_1_variable <= basis_size_2_variable,
        ExcMessage("The second dimension must not be smaller than the first"));

      Assert(quantity == EvaluatorQuantity::value, ExcInternalError());

      // we do recursion until dim==1 or dim==2 and we have
      // basis_size_1==basis_size_2. The latter optimization increases
      // optimization possibilities for the compiler but does only work for
      // aliased pointers if the sizes are equal.
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
                 AlignedVector<Number2>(),
                 AlignedVector<Number2>(),
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
              FEEvaluationImplBasisChange<
                variant,
                quantity,
                next_dim,
                basis_size_1,
                basis_size_2,
                Number,
                Number2>::do_forward(1,
                                     transformation_matrix,
                                     values_in +
                                       (q - 1) *
                                         Utilities::fixed_power<next_dim>(np_1),
                                     values_out +
                                       (q - 1) *
                                         Utilities::fixed_power<next_dim>(np_2),
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
     * @param basis_size_1_variable In case the template argument basis_size_1
     * is zero, the size of the first basis can alternatively be passed in as a
     * run time argument. The template argument takes precedence in case it is
     * nonzero for efficiency reasons.
     * @param basis_size_2_variable In case the template argument basis_size_1
     * is zero, the size of the second basis can alternatively be passed in as a
     * run time argument.
     */
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
    static void
    do_backward(
      const unsigned int            n_components,
      const AlignedVector<Number2> &transformation_matrix,
      const bool                    add_into_result,
      Number *                      values_in,
      Number *                      values_out,
      const unsigned int basis_size_1_variable = numbers::invalid_unsigned_int,
      const unsigned int basis_size_2_variable = numbers::invalid_unsigned_int)
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
                                          basis_size_2,
                                          Number,
                                          Number2>::
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
    static void
    do_mass(const unsigned int            n_components,
            const AlignedVector<Number2> &transformation_matrix,
            const AlignedVector<Number> & coefficients,
            const Number *                values_in,
            Number *                      scratch_data,
            Number *                      values_out)
    {
      constexpr int next_dim = dim > 1 ? dim - 1 : dim;
      Number *      my_scratch =
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
          basis_size_2,
          Number,
          Number2>::do_forward(n_components,
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
        FEEvaluationImplBasisChange<
          variant,
          EvaluatorQuantity::value,
          next_dim,
          basis_size_1,
          basis_size_2,
          Number,
          Number2>::do_backward(n_components,
                                transformation_matrix,
                                false,
                                my_scratch +
                                  q * Utilities::pow(basis_size_2, dim - 1),
                                values_out +
                                  q * Utilities::pow(basis_size_1, dim - 1));
    }
  };



  /**
   * This struct performs the evaluation of function values, gradients and
   * Hessians for tensor-product finite elements. This a specialization for
   * elements where the nodal points coincide with the quadrature points like
   * FE_Q shape functions on Gauss-Lobatto elements integrated with
   * Gauss-Lobatto quadrature. The assumption of this class is that the shape
   * 'values' operation is identity, which allows us to write shorter code.
   *
   * In literature, this form of evaluation is often called spectral
   * evaluation, spectral collocation or simply collocation, meaning the same
   * location for shape functions and evaluation space (quadrature points).
   */
  template <int dim, int fe_degree, typename Number>
  struct FEEvaluationImplCollocation
  {
    static void
    evaluate(const unsigned int                     n_components,
             const EvaluationFlags::EvaluationFlags evaluation_flag,
             const Number *                         values_dofs,
             FEEvaluationData<dim, Number, false> & fe_eval);

    static void
    do_evaluate(const MatrixFreeFunctions::UnivariateShapeData<Number> &shape,
                const EvaluationFlags::EvaluationFlags evaluation_flag,
                const Number *                         values_dofs,
                Number *                               gradients_quad,
                Number *                               hessians_quad);

    static void
    integrate(const unsigned int                     n_components,
              const EvaluationFlags::EvaluationFlags integration_flag,
              Number *                               values_dofs,
              FEEvaluationData<dim, Number, false> & fe_eval,
              const bool                             add_into_values_array);

    static void
    do_integrate(const MatrixFreeFunctions::UnivariateShapeData<Number> &shape,
                 const EvaluationFlags::EvaluationFlags integration_flag,
                 Number *                               values_dofs,
                 Number *                               gradients_quad,
                 const Number *                         hessians_quad,
                 const bool                             add_into_values_array);
  };



  template <int dim, int fe_degree, typename Number>
  inline void
  FEEvaluationImplCollocation<dim, fe_degree, Number>::evaluate(
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags evaluation_flag,
    const Number *                         values_dofs,
    FEEvaluationData<dim, Number, false> & fe_eval)
  {
    constexpr std::size_t n_points = Utilities::pow(fe_degree + 1, dim);

    for (unsigned int c = 0; c < n_components; ++c)
      {
        if ((evaluation_flag & EvaluationFlags::values) != 0u)
          for (unsigned int i = 0; i < n_points; ++i)
            fe_eval.begin_values()[n_points * c + i] =
              values_dofs[n_points * c + i];

        do_evaluate(fe_eval.get_shape_info().data.front(),
                    evaluation_flag,
                    values_dofs + c * n_points,
                    fe_eval.begin_gradients() + c * dim * n_points,
                    fe_eval.begin_hessians() +
                      c * dim * (dim + 1) / 2 * n_points);
      }
  }



  template <int dim, int fe_degree, typename Number>
  inline void
  FEEvaluationImplCollocation<dim, fe_degree, Number>::do_evaluate(
    const MatrixFreeFunctions::UnivariateShapeData<Number> &shape,
    const EvaluationFlags::EvaluationFlags                  evaluation_flag,
    const Number *                                          values_dofs,
    Number *                                                gradients_quad,
    Number *                                                hessians_quad)
  {
    AssertDimension(shape.shape_gradients_collocation_eo.size(),
                    (fe_degree + 2) / 2 * (fe_degree + 1));
    constexpr std::size_t n_points = Utilities::pow(fe_degree + 1, dim);

    EvaluatorTensorProduct<evaluate_evenodd,
                           dim,
                           fe_degree + 1,
                           fe_degree + 1,
                           Number>
      eval(AlignedVector<Number>(),
           shape.shape_gradients_collocation_eo,
           shape.shape_hessians_collocation_eo);
    if ((evaluation_flag &
         (EvaluationFlags::gradients | EvaluationFlags::hessians)) != 0u)
      {
        eval.template gradients<0, true, false>(values_dofs, gradients_quad);
        if (dim > 1)
          eval.template gradients<1, true, false>(values_dofs,
                                                  gradients_quad + n_points);
        if (dim > 2)
          eval.template gradients<2, true, false>(values_dofs,
                                                  gradients_quad +
                                                    2 * n_points);
      }
    if (evaluation_flag & EvaluationFlags::hessians)
      {
        eval.template hessians<0, true, false>(values_dofs, hessians_quad);
        if (dim > 1)
          {
            eval.template gradients<1, true, false>(gradients_quad,
                                                    hessians_quad +
                                                      dim * n_points);
            eval.template hessians<1, true, false>(values_dofs,
                                                   hessians_quad + n_points);
          }
        if (dim > 2)
          {
            eval.template gradients<2, true, false>(gradients_quad,
                                                    hessians_quad +
                                                      4 * n_points);
            eval.template gradients<2, true, false>(gradients_quad + n_points,
                                                    hessians_quad +
                                                      5 * n_points);
            eval.template hessians<2, true, false>(values_dofs,
                                                   hessians_quad +
                                                     2 * n_points);
          }
      }
  }



  template <int dim, int fe_degree, typename Number>
  inline void
  FEEvaluationImplCollocation<dim, fe_degree, Number>::integrate(
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags integration_flag,
    Number *                               values_dofs,
    FEEvaluationData<dim, Number, false> & fe_eval,
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

        do_integrate(fe_eval.get_shape_info().data.front(),
                     integration_flag,
                     values_dofs + c * n_points,
                     fe_eval.begin_gradients() + c * dim * n_points,
                     fe_eval.begin_hessians() +
                       c * dim * (dim + 1) / 2 * n_points,
                     add_into_values_array ||
                       ((integration_flag & EvaluationFlags::values) != 0u));
      }
  }



  template <int dim, int fe_degree, typename Number>
  inline void
  FEEvaluationImplCollocation<dim, fe_degree, Number>::do_integrate(
    const MatrixFreeFunctions::UnivariateShapeData<Number> &shape,
    const EvaluationFlags::EvaluationFlags                  integration_flag,
    Number *                                                values_dofs,
    Number *                                                gradients_quad,
    const Number *                                          hessians_quad,
    const bool add_into_values_array)
  {
    AssertDimension(shape.shape_gradients_collocation_eo.size(),
                    (fe_degree + 2) / 2 * (fe_degree + 1));

    EvaluatorTensorProduct<evaluate_evenodd,
                           dim,
                           fe_degree + 1,
                           fe_degree + 1,
                           Number>
                          eval(AlignedVector<Number>(),
           shape.shape_gradients_collocation_eo,
           shape.shape_hessians_collocation_eo);
    constexpr std::size_t n_points = Utilities::pow(fe_degree + 1, dim);

    if ((integration_flag & EvaluationFlags::hessians) != 0u)
      {
        // diagonal
        // grad xx
        if (add_into_values_array == true)
          eval.template hessians<0, false, true>(hessians_quad, values_dofs);
        else
          eval.template hessians<0, false, false>(hessians_quad, values_dofs);
        // grad yy
        if (dim > 1)
          eval.template hessians<1, false, true>(hessians_quad + n_points,
                                                 values_dofs);
        // grad zz
        if (dim > 2)
          eval.template hessians<2, false, true>(hessians_quad + 2 * n_points,
                                                 values_dofs);
        // off-diagonal
        if (dim == 2)
          {
            // grad xy, queue into gradient
            if (integration_flag & EvaluationFlags::gradients)
              eval.template gradients<1, false, true>(hessians_quad +
                                                        2 * n_points,
                                                      gradients_quad);
            else
              eval.template gradients<1, false, false>(hessians_quad +
                                                         2 * n_points,
                                                       gradients_quad);
          }
        if (dim == 3)
          {
            // grad xy, queue into gradient
            if (integration_flag & EvaluationFlags::gradients)
              eval.template gradients<1, false, true>(hessians_quad +
                                                        3 * n_points,
                                                      gradients_quad);
            else
              eval.template gradients<1, false, false>(hessians_quad +
                                                         3 * n_points,
                                                       gradients_quad);

            // grad xz
            eval.template gradients<2, false, true>(hessians_quad +
                                                      4 * n_points,
                                                    gradients_quad);

            // grad yz
            if (integration_flag & EvaluationFlags::gradients)
              eval.template gradients<2, false, true>(
                hessians_quad + 5 * n_points, gradients_quad + n_points);
            else
              eval.template gradients<2, false, false>(
                hessians_quad + 5 * n_points, gradients_quad + n_points);
          }

        // if we did not integrate gradients, set the last slot to zero
        // which was not touched before, in order to avoid the if
        // statement in the gradients loop below
        if ((integration_flag & EvaluationFlags::gradients) == 0u)
          for (unsigned int q = 0; q < n_points; ++q)
            gradients_quad[(dim - 1) * n_points + q] = Number();
      }

    if ((integration_flag &
         (EvaluationFlags::gradients | EvaluationFlags::hessians)) != 0u)
      {
        if (add_into_values_array ||
            (integration_flag & EvaluationFlags::hessians) != 0u)
          eval.template gradients<0, false, true>(gradients_quad, values_dofs);
        else
          eval.template gradients<0, false, false>(gradients_quad, values_dofs);
        if (dim > 1)
          eval.template gradients<1, false, true>(gradients_quad + n_points,
                                                  values_dofs);
        if (dim > 2)
          eval.template gradients<2, false, true>(gradients_quad + 2 * n_points,
                                                  values_dofs);
      }
  }



  /**
   * This struct performs the evaluation of function values, gradients and
   * Hessians for tensor-product finite elements. This a specialization for
   * symmetric basis functions about the mid point 0.5 of the unit interval
   * with the same number of quadrature points as degrees of freedom. In that
   * case, we can first transform the basis to one that has the nodal points
   * in the quadrature points (i.e., the collocation space) and then perform
   * the evaluation of the first and second derivatives in this transformed
   * space, using the identity operation for the shape values.
   */
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct FEEvaluationImplTransformToCollocation
  {
    static void
    evaluate(const unsigned int                     n_components,
             const EvaluationFlags::EvaluationFlags evaluation_flag,
             const Number *                         values_dofs,
             FEEvaluationData<dim, Number, false> & fe_eval);

    static void
    integrate(const unsigned int                     n_components,
              const EvaluationFlags::EvaluationFlags evaluation_flag,
              Number *                               values_dofs,
              FEEvaluationData<dim, Number, false> & fe_eval,
              const bool                             add_into_values_array);
  };



  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  inline void
  FEEvaluationImplTransformToCollocation<
    dim,
    fe_degree,
    n_q_points_1d,
    Number>::evaluate(const unsigned int                     n_components,
                      const EvaluationFlags::EvaluationFlags evaluation_flag,
                      const Number *                         values_dofs,
                      FEEvaluationData<dim, Number, false> & fe_eval)
  {
    const auto &shape_data = fe_eval.get_shape_info().data.front();

    Assert(n_q_points_1d > fe_degree,
           ExcMessage("You lose information when going to a collocation space "
                      "of lower degree, so the evaluation results would be "
                      "wrong. Thus, this class does not permit the desired "
                      "operation."));
    constexpr std::size_t n_dofs     = Utilities::pow(fe_degree + 1, dim);
    constexpr std::size_t n_q_points = Utilities::pow(n_q_points_1d, dim);

    for (unsigned int c = 0; c < n_components; ++c)
      {
        FEEvaluationImplBasisChange<
          evaluate_evenodd,
          EvaluatorQuantity::value,
          dim,
          (fe_degree >= n_q_points_1d ? n_q_points_1d : fe_degree + 1),
          n_q_points_1d,
          Number,
          Number>::do_forward(1,
                              shape_data.shape_values_eo,
                              values_dofs + c * n_dofs,
                              fe_eval.begin_values() + c * n_q_points);

        // apply derivatives in the collocation space
        if (evaluation_flag &
            (EvaluationFlags::gradients | EvaluationFlags::hessians))
          FEEvaluationImplCollocation<dim, n_q_points_1d - 1, Number>::
            do_evaluate(shape_data,
                        evaluation_flag & (EvaluationFlags::gradients |
                                           EvaluationFlags::hessians),
                        fe_eval.begin_values() + c * n_q_points,
                        fe_eval.begin_gradients() + c * dim * n_q_points,
                        fe_eval.begin_hessians() +
                          c * dim * (dim + 1) / 2 * n_q_points);
      }
  }



  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  inline void
  FEEvaluationImplTransformToCollocation<
    dim,
    fe_degree,
    n_q_points_1d,
    Number>::integrate(const unsigned int                     n_components,
                       const EvaluationFlags::EvaluationFlags integration_flag,
                       Number *                               values_dofs,
                       FEEvaluationData<dim, Number, false> & fe_eval,
                       const bool add_into_values_array)
  {
    const auto &shape_data = fe_eval.get_shape_info().data.front();

    Assert(n_q_points_1d > fe_degree,
           ExcMessage("You lose information when going to a collocation space "
                      "of lower degree, so the evaluation results would be "
                      "wrong. Thus, this class does not permit the desired "
                      "operation."));
    constexpr std::size_t n_q_points = Utilities::pow(n_q_points_1d, dim);

    for (unsigned int c = 0; c < n_components; ++c)
      {
        // apply derivatives in collocation space
        if (integration_flag &
            (EvaluationFlags::gradients | EvaluationFlags::hessians))
          FEEvaluationImplCollocation<dim, n_q_points_1d - 1, Number>::
            do_integrate(shape_data,
                         integration_flag & (EvaluationFlags::gradients |
                                             EvaluationFlags::hessians),
                         fe_eval.begin_values() + c * n_q_points,
                         fe_eval.begin_gradients() + c * dim * n_q_points,
                         fe_eval.begin_hessians() +
                           c * dim * (dim + 1) / 2 * n_q_points,
                         /*add_into_values_array=*/
                         integration_flag & EvaluationFlags::values);

        // transform back to the original space
        FEEvaluationImplBasisChange<
          evaluate_evenodd,
          EvaluatorQuantity::value,
          dim,
          (fe_degree >= n_q_points_1d ? n_q_points_1d : fe_degree + 1),
          n_q_points_1d,
          Number,
          Number>::do_backward(1,
                               shape_data.shape_values_eo,
                               add_into_values_array,
                               fe_eval.begin_values() + c * n_q_points,
                               values_dofs +
                                 c * Utilities::pow(fe_degree + 1, dim));
      }
  }



  /**
   * Helper function to specify whether transformation to collocation should
   * be used: It should give correct results (first condition), we need to be
   * able to initialize the fields in shape_info.templates.h from the
   * polynomials (second condition), and it should be the most efficient
   * choice in terms of operation counts (third condition).
   */
  constexpr bool
  use_collocation_evaluation(const unsigned int fe_degree,
                             const unsigned int n_q_points_1d)
  {
    return (n_q_points_1d > fe_degree) && (n_q_points_1d < 200) &&
           (n_q_points_1d <= 3 * fe_degree / 2 + 1);
  }


  /**
   * This class chooses an appropriate evaluation strategy based on the
   * template parameters and the shape_info variable which contains runtime
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
  template <int dim, typename Number>
  struct FEEvaluationImplEvaluateSelector
  {
    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                     n_components,
        const EvaluationFlags::EvaluationFlags evaluation_flag,
        const Number *                         values_dofs,
        FEEvaluationData<dim, Number, false> & fe_eval)
    {
      const auto element_type = fe_eval.get_shape_info().element_type;
      using ElementType       = MatrixFreeFunctions::ElementType;

      Assert(fe_eval.get_shape_info().data.size() == 1 ||
               (fe_eval.get_shape_info().data.size() == dim &&
                element_type == ElementType::tensor_general) ||
               element_type == ElementType::tensor_raviart_thomas,
             ExcNotImplemented());

      if (fe_degree >= 0 && fe_degree + 1 == n_q_points_1d &&
          element_type == ElementType::tensor_symmetric_collocation)
        {
          FEEvaluationImplCollocation<dim, fe_degree, Number>::evaluate(
            n_components, evaluation_flag, values_dofs, fe_eval);
        }
      // '<=' on type means tensor_symmetric or tensor_symmetric_hermite, see
      // shape_info.h for more details
      else if (fe_degree >= 0 &&
               use_collocation_evaluation(fe_degree, n_q_points_1d) &&
               element_type <= ElementType::tensor_symmetric)
        {
          FEEvaluationImplTransformToCollocation<
            dim,
            fe_degree,
            n_q_points_1d,
            Number>::evaluate(n_components,
                              evaluation_flag,
                              values_dofs,
                              fe_eval);
        }
      else if (fe_degree >= 0 && element_type <= ElementType::tensor_symmetric)
        {
          FEEvaluationImpl<ElementType::tensor_symmetric,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::evaluate(n_components,
                                             evaluation_flag,
                                             values_dofs,
                                             fe_eval);
        }
      else if (element_type == ElementType::tensor_symmetric_plus_dg0)
        {
          FEEvaluationImpl<ElementType::tensor_symmetric_plus_dg0,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::evaluate(n_components,
                                             evaluation_flag,
                                             values_dofs,
                                             fe_eval);
        }
      else if (element_type == ElementType::truncated_tensor)
        {
          FEEvaluationImpl<ElementType::truncated_tensor,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::evaluate(n_components,
                                             evaluation_flag,
                                             values_dofs,
                                             fe_eval);
        }
      else if (element_type == ElementType::tensor_none)
        {
          FEEvaluationImpl<ElementType::tensor_none,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::evaluate(n_components,
                                             evaluation_flag,
                                             values_dofs,
                                             fe_eval);
        }
      else if (element_type == ElementType::tensor_raviart_thomas)
        {
          FEEvaluationImpl<
            ElementType::tensor_raviart_thomas,
            dim,
            (fe_degree == -1) ? 1 : fe_degree,
            (n_q_points_1d < 1) ? 1 : n_q_points_1d,
            Number>::template evaluate_or_integrate<false>(evaluation_flag,
                                                           const_cast<Number *>(
                                                             values_dofs),
                                                           fe_eval);
        }
      else
        {
          FEEvaluationImpl<ElementType::tensor_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::evaluate(n_components,
                                             evaluation_flag,
                                             values_dofs,
                                             fe_eval);
        }

      return false;
    }
  };



  /**
   * This class chooses an appropriate evaluation strategy based on the
   * template parameters and the shape_info variable which contains runtime
   * parameters for the strategy underlying FEEvaluation::integrate(), i.e.
   * this calls internal::FEEvaluationImpl::integrate(),
   * internal::FEEvaluationImplCollocation::integrate() or
   * internal::FEEvaluationImplTransformToCollocation::integrate() with
   * appropriate template parameters. In case the template parameters
   * fe_degree and n_q_points_1d contain valid information (i.e. fe_degree>-1
   * and n_q_points_1d>0), we simply pass these values to the respective
   * template specializations.  Otherwise, we perform a runtime matching of
   * the runtime parameters to find the correct specialization. This matching
   * currently supports $0\leq fe\_degree \leq 9$ and $degree+1\leq
   * n\_q\_points\_1d\leq fe\_degree+2$.
   */
  template <int dim, typename Number>
  struct FEEvaluationImplIntegrateSelector
  {
    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                     n_components,
        const EvaluationFlags::EvaluationFlags integration_flag,
        Number *                               values_dofs,
        FEEvaluationData<dim, Number, false> & fe_eval,
        const bool                             sum_into_values_array)
    {
      const auto element_type = fe_eval.get_shape_info().element_type;
      using ElementType       = MatrixFreeFunctions::ElementType;

      Assert(fe_eval.get_shape_info().data.size() == 1 ||
               (fe_eval.get_shape_info().data.size() == dim &&
                element_type == ElementType::tensor_general) ||
               element_type == ElementType::tensor_raviart_thomas,
             ExcNotImplemented());

      if (fe_degree >= 0 && fe_degree + 1 == n_q_points_1d &&
          element_type == ElementType::tensor_symmetric_collocation)
        {
          FEEvaluationImplCollocation<dim, fe_degree, Number>::integrate(
            n_components,
            integration_flag,
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
          FEEvaluationImplTransformToCollocation<
            dim,
            fe_degree,
            n_q_points_1d,
            Number>::integrate(n_components,
                               integration_flag,
                               values_dofs,
                               fe_eval,
                               sum_into_values_array);
        }
      else if (fe_degree >= 0 && element_type <= ElementType::tensor_symmetric)
        {
          FEEvaluationImpl<ElementType::tensor_symmetric,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::integrate(n_components,
                                              integration_flag,
                                              values_dofs,
                                              fe_eval,
                                              sum_into_values_array);
        }
      else if (element_type == ElementType::tensor_symmetric_plus_dg0)
        {
          FEEvaluationImpl<ElementType::tensor_symmetric_plus_dg0,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::integrate(n_components,
                                              integration_flag,
                                              values_dofs,
                                              fe_eval,
                                              sum_into_values_array);
        }
      else if (element_type == ElementType::truncated_tensor)
        {
          FEEvaluationImpl<ElementType::truncated_tensor,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::integrate(n_components,
                                              integration_flag,
                                              values_dofs,
                                              fe_eval,
                                              sum_into_values_array);
        }
      else if (element_type == ElementType::tensor_none)
        {
          FEEvaluationImpl<ElementType::tensor_none,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::integrate(n_components,
                                              integration_flag,
                                              values_dofs,
                                              fe_eval,
                                              sum_into_values_array);
        }
      else if (element_type == ElementType::tensor_raviart_thomas)
        {
          FEEvaluationImpl<ElementType::tensor_raviart_thomas,
                           dim,
                           (fe_degree == -1) ? 1 : fe_degree,
                           (n_q_points_1d < 1) ? 1 : n_q_points_1d,
                           Number>::
            template evaluate_or_integrate<true>(integration_flag,
                                                 const_cast<Number *>(
                                                   values_dofs),
                                                 fe_eval,
                                                 sum_into_values_array);
        }
      else
        {
          FEEvaluationImpl<ElementType::tensor_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::integrate(n_components,
                                              integration_flag,
                                              values_dofs,
                                              fe_eval,
                                              sum_into_values_array);
        }

      return false;
    }
  };



  template <bool symmetric_evaluate,
            int  dim,
            int  fe_degree,
            int  n_q_points_1d,
            typename Number>
  struct FEFaceEvaluationImpl
  {
    // We enable a transformation to collocation for derivatives if it gives
    // correct results (first two conditions), if it is the most efficient
    // choice in terms of operation counts (third condition) and if we were
    // able to initialize the fields in shape_info.templates.h from the
    // polynomials (fourth condition).
    using Eval = EvaluatorTensorProduct<symmetric_evaluate ? evaluate_evenodd :
                                                             evaluate_general,
                                        dim - 1,
                                        fe_degree + 1,
                                        n_q_points_1d,
                                        Number>;

    static Eval
    create_evaluator_tensor_product(
      const MatrixFreeFunctions::UnivariateShapeData<Number> &data,
      const unsigned int                                      subface_index,
      const unsigned int                                      direction)
    {
      if (symmetric_evaluate)
        return Eval(data.shape_values_eo,
                    data.shape_gradients_eo,
                    data.shape_hessians_eo,
                    data.fe_degree + 1,
                    data.n_q_points_1d);
      else if (subface_index >= GeometryInfo<dim>::max_children_per_cell)
        return Eval(data.shape_values,
                    data.shape_gradients,
                    data.shape_hessians,
                    data.fe_degree + 1,
                    data.n_q_points_1d);
      else
        {
          const unsigned int index =
            direction == 0 ? subface_index % 2 : subface_index / 2;
          return Eval(data.values_within_subface[index],
                      data.gradients_within_subface[index],
                      data.hessians_within_subface[index],
                      data.fe_degree + 1,
                      data.n_q_points_1d);
        }
    }

    static void
    evaluate_in_face(
      const unsigned int                                      n_components,
      const EvaluationFlags::EvaluationFlags                  evaluation_flag,
      const MatrixFreeFunctions::UnivariateShapeData<Number> &data,
      Number *                                                values_dofs,
      Number *                                                values_quad,
      Number *                                                gradients_quad,
      Number *                                                hessians_quad,
      Number *                                                scratch_data,
      const unsigned int                                      subface_index)
    {
      Eval eval0 = create_evaluator_tensor_product(data, subface_index, 0);
      Eval eval1 = create_evaluator_tensor_product(data, subface_index, 1);

      const std::size_t n_dofs = fe_degree > -1 ?
                                   Utilities::pow(fe_degree + 1, dim - 1) :
                                   Utilities::pow(data.fe_degree + 1, dim - 1);
      const std::size_t n_q_points =
        fe_degree > -1 ? Utilities::pow(n_q_points_1d, dim - 1) :
                         Utilities::pow(data.n_q_points_1d, dim - 1);

      // keep a copy of the original pointer for the case of the Hessians
      Number *values_dofs_ptr = values_dofs;

      if ((evaluation_flag & EvaluationFlags::values) != 0u &&
          ((evaluation_flag & EvaluationFlags::gradients) == 0u))
        for (unsigned int c = 0; c < n_components; ++c)
          {
            switch (dim)
              {
                case 3:
                  eval0.template values<0, true, false>(values_dofs,
                                                        values_quad);
                  eval1.template values<1, true, false>(values_quad,
                                                        values_quad);
                  break;
                case 2:
                  eval0.template values<0, true, false>(values_dofs,
                                                        values_quad);
                  break;
                case 1:
                  values_quad[0] = values_dofs[0];
                  break;
                default:
                  Assert(false, ExcNotImplemented());
              }
            // Note: we always keep storage of values, 1st and 2nd derivatives
            // in an array
            values_dofs += 3 * n_dofs;
            values_quad += n_q_points;
          }
      else if ((evaluation_flag & EvaluationFlags::gradients) != 0u)
        for (unsigned int c = 0; c < n_components; ++c)
          {
            switch (dim)
              {
                case 3:
                  if (symmetric_evaluate &&
                      use_collocation_evaluation(fe_degree, n_q_points_1d))
                    {
                      eval0.template values<0, true, false>(values_dofs,
                                                            values_quad);
                      eval0.template values<1, true, false>(values_quad,
                                                            values_quad);
                      EvaluatorTensorProduct<evaluate_evenodd,
                                             dim - 1,
                                             n_q_points_1d,
                                             n_q_points_1d,
                                             Number>
                        eval_grad(AlignedVector<Number>(),
                                  data.shape_gradients_collocation_eo,
                                  AlignedVector<Number>());
                      eval_grad.template gradients<0, true, false>(
                        values_quad, gradients_quad);
                      eval_grad.template gradients<1, true, false>(
                        values_quad, gradients_quad + n_q_points);
                    }
                  else
                    {
                      // grad x
                      eval0.template gradients<0, true, false>(values_dofs,
                                                               scratch_data);
                      eval1.template values<1, true, false>(scratch_data,
                                                            gradients_quad);

                      // grad y
                      eval0.template values<0, true, false>(values_dofs,
                                                            scratch_data);
                      eval1.template gradients<1, true, false>(scratch_data,
                                                               gradients_quad +
                                                                 n_q_points);

                      if ((evaluation_flag & EvaluationFlags::values) != 0u)
                        eval1.template values<1, true, false>(scratch_data,
                                                              values_quad);
                    }
                  // grad z
                  eval0.template values<0, true, false>(values_dofs + n_dofs,
                                                        scratch_data);
                  eval1.template values<1, true, false>(
                    scratch_data, gradients_quad + (dim - 1) * n_q_points);

                  break;
                case 2:
                  eval0.template values<0, true, false>(values_dofs + n_dofs,
                                                        gradients_quad +
                                                          n_q_points);
                  eval0.template gradients<0, true, false>(values_dofs,
                                                           gradients_quad);
                  if ((evaluation_flag & EvaluationFlags::values) != 0u)
                    eval0.template values<0, true, false>(values_dofs,
                                                          values_quad);
                  break;
                case 1:
                  values_quad[0]    = values_dofs[0];
                  gradients_quad[0] = values_dofs[1];
                  break;
                default:
                  AssertThrow(false, ExcNotImplemented());
              }
            values_dofs += 3 * n_dofs;
            values_quad += n_q_points;
            gradients_quad += dim * n_q_points;
          }

      if ((evaluation_flag & EvaluationFlags::hessians) != 0u)
        {
          values_dofs = values_dofs_ptr;
          for (unsigned int c = 0; c < n_components; ++c)
            {
              switch (dim)
                {
                  case 3:
                    // grad xx
                    eval0.template hessians<0, true, false>(values_dofs,
                                                            scratch_data);
                    eval1.template values<1, true, false>(scratch_data,
                                                          hessians_quad);

                    // grad yy
                    eval0.template values<0, true, false>(values_dofs,
                                                          scratch_data);
                    eval1.template hessians<1, true, false>(scratch_data,
                                                            hessians_quad +
                                                              n_q_points);

                    // grad zz
                    eval0.template values<0, true, false>(values_dofs +
                                                            2 * n_dofs,
                                                          scratch_data);
                    eval1.template values<1, true, false>(scratch_data,
                                                          hessians_quad +
                                                            2 * n_q_points);

                    // grad xy
                    eval0.template gradients<0, true, false>(values_dofs,
                                                             scratch_data);
                    eval1.template gradients<1, true, false>(scratch_data,
                                                             hessians_quad +
                                                               3 * n_q_points);

                    // grad xz
                    eval0.template gradients<0, true, false>(values_dofs +
                                                               n_dofs,
                                                             scratch_data);
                    eval1.template values<1, true, false>(scratch_data,
                                                          hessians_quad +
                                                            4 * n_q_points);

                    // grad yz
                    eval0.template values<0, true, false>(values_dofs + n_dofs,
                                                          scratch_data);
                    eval1.template gradients<1, true, false>(scratch_data,
                                                             hessians_quad +
                                                               5 * n_q_points);

                    break;
                  case 2:
                    // grad xx
                    eval0.template hessians<0, true, false>(values_dofs,
                                                            hessians_quad);
                    // grad yy
                    eval0.template values<0, true, false>(
                      values_dofs + 2 * n_dofs, hessians_quad + n_q_points);
                    // grad xy
                    eval0.template gradients<0, true, false>(
                      values_dofs + n_dofs, hessians_quad + 2 * n_q_points);
                    break;
                  case 1:
                    hessians_quad[0] = values_dofs[2];
                    break;
                  default:
                    AssertThrow(false, ExcNotImplemented());
                }
              values_dofs += 3 * n_dofs;
              hessians_quad += dim * (dim + 1) / 2 * n_q_points;
            }
        }
    }

    static void
    integrate_in_face(
      const unsigned int                                      n_components,
      const EvaluationFlags::EvaluationFlags                  integration_flag,
      const MatrixFreeFunctions::UnivariateShapeData<Number> &data,
      Number *                                                values_dofs,
      Number *                                                values_quad,
      Number *                                                gradients_quad,
      Number *                                                hessians_quad,
      Number *                                                scratch_data,
      const unsigned int                                      subface_index)
    {
      Eval eval0 = create_evaluator_tensor_product(data, subface_index, 0);
      Eval eval1 = create_evaluator_tensor_product(data, subface_index, 1);

      const std::size_t n_dofs =
        fe_degree > -1 ?
          Utilities::pow(fe_degree + 1, dim - 1) :
          (dim > 1 ? Utilities::fixed_power<dim - 1>(data.fe_degree + 1) : 1);
      const std::size_t n_q_points =
        fe_degree > -1 ? Utilities::pow(n_q_points_1d, dim - 1) :
                         Utilities::pow(data.n_q_points_1d, dim - 1);

      // keep a copy of the original pointer for the case of the Hessians
      Number *values_dofs_ptr = values_dofs;

      if ((integration_flag & EvaluationFlags::values) != 0u &&
          (integration_flag & EvaluationFlags::gradients) == 0u)
        for (unsigned int c = 0; c < n_components; ++c)
          {
            switch (dim)
              {
                case 3:
                  eval1.template values<1, false, false>(values_quad,
                                                         values_quad);
                  eval0.template values<0, false, false>(values_quad,
                                                         values_dofs);
                  break;
                case 2:
                  eval0.template values<0, false, false>(values_quad,
                                                         values_dofs);
                  break;
                case 1:
                  values_dofs[0] = values_quad[0];
                  break;
                default:
                  Assert(false, ExcNotImplemented());
              }
            values_dofs += 3 * n_dofs;
            values_quad += n_q_points;
          }
      else if ((integration_flag & EvaluationFlags::gradients) != 0u)
        for (unsigned int c = 0; c < n_components; ++c)
          {
            switch (dim)
              {
                case 3:
                  // grad z
                  eval1.template values<1, false, false>(gradients_quad +
                                                           2 * n_q_points,
                                                         gradients_quad +
                                                           2 * n_q_points);
                  eval0.template values<0, false, false>(gradients_quad +
                                                           2 * n_q_points,
                                                         values_dofs + n_dofs);
                  if (symmetric_evaluate &&
                      use_collocation_evaluation(fe_degree, n_q_points_1d))
                    {
                      EvaluatorTensorProduct<evaluate_evenodd,
                                             dim - 1,
                                             n_q_points_1d,
                                             n_q_points_1d,
                                             Number>
                        eval_grad(AlignedVector<Number>(),
                                  data.shape_gradients_collocation_eo,
                                  AlignedVector<Number>());
                      if ((integration_flag & EvaluationFlags::values) != 0u)
                        eval_grad.template gradients<1, false, true>(
                          gradients_quad + n_q_points, values_quad);
                      else
                        eval_grad.template gradients<1, false, false>(
                          gradients_quad + n_q_points, values_quad);
                      eval_grad.template gradients<0, false, true>(
                        gradients_quad, values_quad);
                      eval0.template values<1, false, false>(values_quad,
                                                             values_quad);
                      eval0.template values<0, false, false>(values_quad,
                                                             values_dofs);
                    }
                  else
                    {
                      if ((integration_flag & EvaluationFlags::values) != 0u)
                        {
                          eval1.template values<1, false, false>(values_quad,
                                                                 scratch_data);
                          eval1.template gradients<1, false, true>(
                            gradients_quad + n_q_points, scratch_data);
                        }
                      else
                        eval1.template gradients<1, false, false>(
                          gradients_quad + n_q_points, scratch_data);

                      // grad y
                      eval0.template values<0, false, false>(scratch_data,
                                                             values_dofs);

                      // grad x
                      eval1.template values<1, false, false>(gradients_quad,
                                                             scratch_data);
                      eval0.template gradients<0, false, true>(scratch_data,
                                                               values_dofs);
                    }
                  break;
                case 2:
                  eval0.template values<0, false, false>(gradients_quad +
                                                           n_q_points,
                                                         values_dofs + n_dofs);
                  eval0.template gradients<0, false, false>(gradients_quad,
                                                            values_dofs);
                  if ((integration_flag & EvaluationFlags::values) != 0u)
                    eval0.template values<0, false, true>(values_quad,
                                                          values_dofs);
                  break;
                case 1:
                  values_dofs[0] = values_quad[0];
                  values_dofs[1] = gradients_quad[0];
                  break;
                default:
                  AssertThrow(false, ExcNotImplemented());
              }
            values_dofs += 3 * n_dofs;
            values_quad += n_q_points;
            gradients_quad += dim * n_q_points;
          }

      if ((integration_flag & EvaluationFlags::hessians) != 0u)
        {
          values_dofs = values_dofs_ptr;
          for (unsigned int c = 0; c < n_components; ++c)
            {
              switch (dim)
                {
                  case 3:
                    // grad xx
                    eval1.template values<1, false, false>(hessians_quad,
                                                           scratch_data);
                    if ((integration_flag & (EvaluationFlags::values |
                                             EvaluationFlags::gradients)) != 0u)
                      eval0.template hessians<0, false, true>(scratch_data,
                                                              values_dofs);
                    else
                      eval0.template hessians<0, false, false>(scratch_data,
                                                               values_dofs);

                    // grad yy
                    eval1.template hessians<1, false, false>(hessians_quad +
                                                               n_q_points,
                                                             scratch_data);
                    eval0.template values<0, false, true>(scratch_data,
                                                          values_dofs);

                    // grad zz
                    eval1.template values<1, false, false>(hessians_quad +
                                                             2 * n_q_points,
                                                           scratch_data);
                    eval0.template values<0, false, false>(scratch_data,
                                                           values_dofs +
                                                             2 * n_dofs);

                    // grad xy
                    eval1.template gradients<1, false, false>(hessians_quad +
                                                                3 * n_q_points,
                                                              scratch_data);
                    eval0.template gradients<0, false, true>(scratch_data,
                                                             values_dofs);

                    // grad xz
                    eval1.template values<1, false, false>(hessians_quad +
                                                             4 * n_q_points,
                                                           scratch_data);
                    if ((integration_flag & EvaluationFlags::gradients) != 0u)
                      eval0.template gradients<0, false, true>(scratch_data,
                                                               values_dofs +
                                                                 n_dofs);
                    else
                      eval0.template gradients<0, false, false>(scratch_data,
                                                                values_dofs +
                                                                  n_dofs);

                    // grad yz
                    eval1.template gradients<1, false, false>(hessians_quad +
                                                                5 * n_q_points,
                                                              scratch_data);
                    eval0.template values<0, false, true>(scratch_data,
                                                          values_dofs + n_dofs);

                    break;
                  case 2:
                    // grad xx
                    if ((integration_flag & (EvaluationFlags::values |
                                             EvaluationFlags::gradients)) != 0u)
                      eval0.template hessians<0, false, true>(hessians_quad,
                                                              values_dofs);
                    else
                      eval0.template hessians<0, false, false>(hessians_quad,
                                                               values_dofs);

                    // grad yy
                    eval0.template values<0, false, false>(
                      hessians_quad + n_q_points, values_dofs + 2 * n_dofs);
                    // grad xy
                    if ((integration_flag & EvaluationFlags::gradients) != 0u)
                      eval0.template gradients<0, false, true>(
                        hessians_quad + 2 * n_q_points, values_dofs + n_dofs);
                    else
                      eval0.template gradients<0, false, false>(
                        hessians_quad + 2 * n_q_points, values_dofs + n_dofs);
                    break;
                  case 1:
                    values_dofs[2] = hessians_quad[0];
                    if ((integration_flag & EvaluationFlags::values) == 0u)
                      values_dofs[0] = 0;
                    if ((integration_flag & EvaluationFlags::gradients) == 0u)
                      values_dofs[1] = 0;
                    break;
                  default:
                    AssertThrow(false, ExcNotImplemented());
                }
              values_dofs += 3 * n_dofs;
              hessians_quad += dim * (dim + 1) / 2 * n_q_points;
            }
        }
    }
  };

  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct FEFaceEvaluationImplRaviartThomas
  {
    using EvalGeneral = EvaluatorTensorProduct<evaluate_general,
                                               dim - 1,
                                               fe_degree,
                                               n_q_points_1d,
                                               Number>;
    template <typename EvalType>
    static EvalType
    create_evaluator_tensor_product(
      const MatrixFreeFunctions::UnivariateShapeData<Number> &data,
      const unsigned int                                      subface_index,
      const unsigned int                                      direction)
    {
      if (subface_index >= GeometryInfo<dim>::max_children_per_cell)
        return EvalType(data.shape_values,
                        data.shape_gradients,
                        data.shape_hessians);
      else
        {
          const unsigned int index =
            direction == 0 ? subface_index % 2 : subface_index / 2;
          return EvalType(data.values_within_subface[index],
                          data.gradients_within_subface[index],
                          data.hessians_within_subface[index]);
        }
    }

    template <bool integrate>
    static void
    evaluate_or_integrate_in_face(
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      Number *                               values_dofs,
      FEEvaluationData<dim, Number, true> &  fe_eval,
      Number *                               scratch_data,
      const unsigned int                     subface_index,
      const unsigned int                     face_no)
    {
      const unsigned int face_direction = face_no / 2;

      // We first evaluate the anisotropic faces, i.e the faces where
      // face_direction != component. Note that the call order here is not
      // important, since the pointers are shifted accordingly within the
      // function. However, this is the order in which the components will be in
      // the quadrature points. Furthermore, the isotropic faces have no "normal
      // direction" but we still pass in normal_dir = 2 since this is used for
      // the pointers.
      // -----------------------------------------------------------------------------------
      // |          |                   Anisotropic faces                 | Isotropic faces|
      // | Face dir | comp, coords, normal_dir | comp, coords, normal_dir | comp, coords   |
      // | --------------------------------------------------------------------------------|
      // |    0     | 1, y, 0                  | -                        | 0, y           |
      // |    1     | 0, x, 0                  | -                        | 1, x           |
      // | --------------------------------------------------------------------------------|
      // |    0     | 1, yz, 0                 | 2, yz, 1                 | 0, yz          |
      // |    1     | 2, zx, 0                 | 0, zx, 1                 | 1, zx          |
      // |    2     | 0, xy, 0                 | 1, xy, 1                 | 2, xy          |
      // -----------------------------------------------------------------------------------
      evaluate_in_face_apply<0>(values_dofs,
                                fe_eval,
                                scratch_data,
                                evaluation_flag,
                                face_direction,
                                subface_index,
                                std::integral_constant<bool, integrate>());

      if (dim == 3)
        evaluate_in_face_apply<1>(values_dofs,
                                  fe_eval,
                                  scratch_data,
                                  evaluation_flag,
                                  face_direction,
                                  subface_index,
                                  std::integral_constant<bool, integrate>());

      evaluate_in_face_apply<2>(values_dofs,
                                fe_eval,
                                scratch_data,
                                evaluation_flag,
                                face_direction,
                                subface_index,
                                std::integral_constant<bool, integrate>());
    }

    /*
     * Helper function which applies the 1D kernels for on one
     * component in a face. normal_dir indicates the direction of the continuous
     * component of the RT space. std::integral_constant<bool, false> is the
     * evaluation path, and std::integral_constant<bool, true> below is the
     * integration path. These two functions can be fused together since all
     * offsets and pointers are the exact same.
     */
    template <int normal_dir>
    static inline void
    evaluate_in_face_apply(
      Number *                               values_dofs,
      FEEvaluationData<dim, Number, true> &  fe_eval,
      Number *                               scratch_data,
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      const unsigned int                     face_direction,
      const unsigned int                     subface_index,
      std::integral_constant<bool, false>)
    {
      using EvalNormal =
        EvaluatorTensorProductAnisotropic<evaluate_raviart_thomas,
                                          dim - 1,
                                          (fe_degree == -1) ? 1 : fe_degree + 1,
                                          n_q_points_1d,
                                          Number,
                                          normal_dir>;
      using EvalTangent =
        EvaluatorTensorProductAnisotropic<evaluate_raviart_thomas,
                                          dim - 1,
                                          (fe_degree == -1) ? 1 : fe_degree,
                                          n_q_points_1d,
                                          Number,
                                          normal_dir>;

      using TempEval0 = typename std::
        conditional<normal_dir == 0, EvalNormal, EvalTangent>::type;
      using TempEval1 = typename std::
        conditional<normal_dir == 0, EvalTangent, EvalNormal>::type;
      using Eval0 = typename std::
        conditional<normal_dir == 2, EvalGeneral, TempEval0>::type;
      using Eval1 = typename std::
        conditional<normal_dir == 2, EvalGeneral, TempEval1>::type;

      const MatrixFreeFunctions::ShapeInfo<Number> &shape_info =
        fe_eval.get_shape_info();
      Eval0 eval0 = create_evaluator_tensor_product<Eval0>(
        ((normal_dir == 0) ? shape_info.data[0] : shape_info.data[1]),
        subface_index,
        0);
      Eval1 eval1 = create_evaluator_tensor_product<Eval1>(
        ((normal_dir == 1) ? shape_info.data[0] : shape_info.data[1]),
        subface_index,
        1);

      constexpr std::size_t n_q_points = Utilities::pow(n_q_points_1d, dim - 1);
      const std::size_t n_dofs_tangent = shape_info.dofs_per_component_on_face;
      const std::size_t n_dofs_normal =
        n_dofs_tangent - Utilities::pow(fe_degree, dim - 2);
      const std::size_t dofs_stride =
        (std::is_same<Eval0, EvalGeneral>::value) ? n_dofs_normal :
                                                    n_dofs_tangent;

      static constexpr dealii::ndarray<unsigned int, 3, 3> component_table = {
        {{{1, 2, 0}}, {{2, 0, 1}}, {{0, 1, 2}}}};
      const unsigned int component =
        (dim == 2 && normal_dir == 0 && face_direction == 1) ?
          0 :
          component_table[face_direction][normal_dir];

      // Initial offsets
      values_dofs +=
        3 * ((component == 0) ?
               0 :
               ((component == 1) ?
                  ((face_direction == 0) ? n_dofs_normal : n_dofs_tangent) :
                  ((face_direction == 2) ? n_dofs_tangent + n_dofs_tangent :
                                           n_dofs_normal + n_dofs_tangent)));
      const unsigned int shift = (dim == 2) ? normal_dir / 2 : normal_dir;
      Number *values_quad      = fe_eval.begin_values() + n_q_points * shift;
      Number *gradients_quad =
        fe_eval.begin_gradients() + dim * n_q_points * shift;
      Number *hessians_quad =
        fe_eval.begin_hessians() + dim * (dim + 1) / 2 * n_q_points * shift;

      // Evaluation path
      if ((evaluation_flag & EvaluationFlags::values) &&
          !(evaluation_flag & EvaluationFlags::gradients))
        {
          switch (dim)
            {
              case 3:
                eval0.template values<0, true, false>(values_dofs, values_quad);
                eval1.template values<1, true, false>(values_quad, values_quad);
                break;
              case 2:
                eval0.template values<0, true, false>(values_dofs, values_quad);
                break;
              default:
                Assert(false, ExcNotImplemented());
            }
        }
      else if (evaluation_flag & EvaluationFlags::gradients)
        {
          switch (dim)
            {
              case 3:
                // grad x
                eval0.template gradients<0, true, false>(values_dofs,
                                                         scratch_data);
                eval1.template values<1, true, false>(scratch_data,
                                                      gradients_quad);

                // grad y
                eval0.template values<0, true, false>(values_dofs,
                                                      scratch_data);
                eval1.template gradients<1, true, false>(scratch_data,
                                                         gradients_quad +
                                                           n_q_points);

                if (evaluation_flag & EvaluationFlags::values)
                  eval1.template values<1, true, false>(scratch_data,
                                                        values_quad);

                // grad z
                eval0.template values<0, true, false>(values_dofs + dofs_stride,
                                                      scratch_data);
                eval1.template values<1, true, false>(scratch_data,
                                                      gradients_quad +
                                                        2 * n_q_points);

                break;
              case 2:
                eval0.template values<0, true, false>(values_dofs + dofs_stride,
                                                      gradients_quad +
                                                        n_q_points);
                eval0.template gradients<0, true, false>(values_dofs,
                                                         gradients_quad);
                if ((evaluation_flag & EvaluationFlags::values))
                  eval0.template values<0, true, false>(values_dofs,
                                                        values_quad);
                break;
              default:
                AssertThrow(false, ExcNotImplemented());
            }
        }

      if (evaluation_flag & EvaluationFlags::hessians)
        {
          switch (dim)
            {
              case 3:
                // grad xx
                eval0.template hessians<0, true, false>(values_dofs,
                                                        scratch_data);
                eval1.template values<1, true, false>(scratch_data,
                                                      hessians_quad);

                // grad yy
                eval0.template values<0, true, false>(values_dofs,
                                                      scratch_data);
                eval1.template hessians<1, true, false>(scratch_data,
                                                        hessians_quad +
                                                          n_q_points);

                // grad zz
                eval0.template values<0, true, false>(values_dofs +
                                                        2 * dofs_stride,
                                                      scratch_data);
                eval1.template values<1, true, false>(scratch_data,
                                                      hessians_quad +
                                                        2 * n_q_points);

                // grad xy
                eval0.template gradients<0, true, false>(values_dofs,
                                                         scratch_data);
                eval1.template gradients<1, true, false>(scratch_data,
                                                         hessians_quad +
                                                           3 * n_q_points);

                // grad xz
                eval0.template gradients<0, true, false>(values_dofs +
                                                           dofs_stride,
                                                         scratch_data);
                eval1.template values<1, true, false>(scratch_data,
                                                      hessians_quad +
                                                        4 * n_q_points);

                // grad yz
                eval0.template values<0, true, false>(values_dofs + dofs_stride,
                                                      scratch_data);
                eval1.template gradients<1, true, false>(scratch_data,
                                                         hessians_quad +
                                                           5 * n_q_points);

                break;
              case 2:
                // grad xx
                eval0.template hessians<0, true, false>(values_dofs,
                                                        hessians_quad);
                // grad yy
                eval0.template values<0, true, false>(
                  values_dofs + 2 * dofs_stride, hessians_quad + n_q_points);
                // grad xy
                eval0.template gradients<0, true, false>(
                  values_dofs + dofs_stride, hessians_quad + 2 * n_q_points);
                break;
              default:
                AssertThrow(false, ExcNotImplemented());
            }
        }
    }

    template <int normal_dir>
    static inline void
    evaluate_in_face_apply(
      Number *                               values_dofs,
      FEEvaluationData<dim, Number, true> &  fe_eval,
      Number *                               scratch_data,
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      const unsigned int                     face_direction,
      const unsigned int                     subface_index,
      std::integral_constant<bool, true>)
    {
      using EvalNormal =
        EvaluatorTensorProductAnisotropic<evaluate_raviart_thomas,
                                          dim - 1,
                                          (fe_degree == -1) ? 1 : fe_degree + 1,
                                          n_q_points_1d,
                                          Number,
                                          normal_dir>;
      using EvalTangent =
        EvaluatorTensorProductAnisotropic<evaluate_raviart_thomas,
                                          dim - 1,
                                          (fe_degree == -1) ? 1 : fe_degree,
                                          n_q_points_1d,
                                          Number,
                                          normal_dir>;

      using TempEval0 = typename std::
        conditional<normal_dir == 0, EvalNormal, EvalTangent>::type;
      using TempEval1 = typename std::
        conditional<normal_dir == 0, EvalTangent, EvalNormal>::type;
      using Eval0 = typename std::
        conditional<normal_dir == 2, EvalGeneral, TempEval0>::type;
      using Eval1 = typename std::
        conditional<normal_dir == 2, EvalGeneral, TempEval1>::type;

      const MatrixFreeFunctions::ShapeInfo<Number> &shape_info =
        fe_eval.get_shape_info();
      Eval0 eval0 = create_evaluator_tensor_product<Eval0>(
        ((normal_dir == 0) ? shape_info.data[0] : shape_info.data[1]),
        subface_index,
        0);
      Eval1 eval1 = create_evaluator_tensor_product<Eval1>(
        ((normal_dir == 1) ? shape_info.data[0] : shape_info.data[1]),
        subface_index,
        1);

      constexpr std::size_t n_q_points = Utilities::pow(n_q_points_1d, dim - 1);
      const std::size_t n_dofs_tangent = shape_info.dofs_per_component_on_face;
      const std::size_t n_dofs_normal =
        n_dofs_tangent - Utilities::pow(fe_degree, dim - 2);
      const std::size_t dofs_stride =
        (std::is_same<Eval0, EvalGeneral>::value) ? n_dofs_normal :
                                                    n_dofs_tangent;

      static constexpr dealii::ndarray<unsigned int, 3, 3> component_table = {
        {{{1, 2, 0}}, {{2, 0, 1}}, {{0, 1, 2}}}};
      const unsigned int component =
        (dim == 2 && normal_dir == 0 && face_direction == 1) ?
          0 :
          component_table[face_direction][normal_dir];

      // Initial offsets
      values_dofs +=
        3 * ((component == 0) ?
               0 :
               ((component == 1) ?
                  ((face_direction == 0) ? n_dofs_normal : n_dofs_tangent) :
                  ((face_direction == 2) ? n_dofs_tangent + n_dofs_tangent :
                                           n_dofs_normal + n_dofs_tangent)));
      const unsigned int shift = (dim == 2) ? normal_dir / 2 : normal_dir;
      Number *values_quad      = fe_eval.begin_values() + n_q_points * shift;
      Number *gradients_quad =
        fe_eval.begin_gradients() + dim * n_q_points * shift;
      Number *hessians_quad =
        fe_eval.begin_hessians() + dim * (dim + 1) / 2 * n_q_points * shift;

      // Integration path
      if ((evaluation_flag & EvaluationFlags::values) &&
          !(evaluation_flag & EvaluationFlags::gradients))
        {
          switch (dim)
            {
              case 3:
                eval1.template values<1, false, false>(values_quad,
                                                       values_quad);
                eval0.template values<0, false, false>(values_quad,
                                                       values_dofs);
                break;
              case 2:
                eval0.template values<0, false, false>(values_quad,
                                                       values_dofs);
                break;
              default:
                Assert(false, ExcNotImplemented());
            }
        }
      else if (evaluation_flag & EvaluationFlags::gradients)
        {
          switch (dim)
            {
              case 3:
                // grad z
                eval1.template values<1, false, false>(gradients_quad +
                                                         2 * n_q_points,
                                                       gradients_quad +
                                                         2 * n_q_points);
                eval0.template values<0, false, false>(
                  gradients_quad + 2 * n_q_points, values_dofs + dofs_stride);

                if (evaluation_flag & EvaluationFlags::values)
                  {
                    eval1.template values<1, false, false>(values_quad,
                                                           scratch_data);
                    eval1.template gradients<1, false, true>(gradients_quad +
                                                               n_q_points,
                                                             scratch_data);
                  }
                else
                  eval1.template gradients<1, false, false>(gradients_quad +
                                                              n_q_points,
                                                            scratch_data);

                // grad y
                eval0.template values<0, false, false>(scratch_data,
                                                       values_dofs);

                // grad x
                eval1.template values<1, false, false>(gradients_quad,
                                                       scratch_data);
                eval0.template gradients<0, false, true>(scratch_data,
                                                         values_dofs);

                break;
              case 2:
                eval0.template values<0, false, false>(
                  gradients_quad + n_q_points, values_dofs + dofs_stride);
                eval0.template gradients<0, false, false>(gradients_quad,
                                                          values_dofs);
                if (evaluation_flag & EvaluationFlags::values)
                  eval0.template values<0, false, true>(values_quad,
                                                        values_dofs);
                break;
              default:
                AssertThrow(false, ExcNotImplemented());
            }
        }

      if (evaluation_flag & EvaluationFlags::hessians)
        {
          switch (dim)
            {
              case 3:
                // grad xx
                eval1.template values<1, false, false>(hessians_quad,
                                                       scratch_data);
                if ((evaluation_flag &
                     (EvaluationFlags::values | EvaluationFlags::gradients)))
                  eval0.template hessians<0, false, true>(scratch_data,
                                                          values_dofs);
                else
                  eval0.template hessians<0, false, false>(scratch_data,
                                                           values_dofs);

                // grad yy
                eval1.template hessians<1, false, false>(hessians_quad +
                                                           n_q_points,
                                                         scratch_data);
                eval0.template values<0, false, true>(scratch_data,
                                                      values_dofs);

                // grad zz
                eval1.template values<1, false, false>(hessians_quad +
                                                         2 * n_q_points,
                                                       scratch_data);
                eval0.template values<0, false, false>(scratch_data,
                                                       values_dofs +
                                                         2 * dofs_stride);

                // grad xy
                eval1.template gradients<1, false, false>(hessians_quad +
                                                            3 * n_q_points,
                                                          scratch_data);
                eval0.template gradients<0, false, true>(scratch_data,
                                                         values_dofs);

                // grad xz
                eval1.template values<1, false, false>(hessians_quad +
                                                         4 * n_q_points,
                                                       scratch_data);
                if ((evaluation_flag & EvaluationFlags::gradients))
                  eval0.template gradients<0, false, true>(scratch_data,
                                                           values_dofs +
                                                             dofs_stride);
                else
                  eval0.template gradients<0, false, false>(scratch_data,
                                                            values_dofs +
                                                              dofs_stride);

                // grad yz
                eval1.template gradients<1, false, false>(hessians_quad +
                                                            5 * n_q_points,
                                                          scratch_data);
                eval0.template values<0, false, true>(scratch_data,
                                                      values_dofs +
                                                        dofs_stride);

                break;
              case 2:
                // grad xx
                if (evaluation_flag &
                    (EvaluationFlags::values | EvaluationFlags::gradients))
                  eval0.template hessians<0, false, true>(hessians_quad,
                                                          values_dofs);
                else
                  eval0.template hessians<0, false, false>(hessians_quad,
                                                           values_dofs);

                // grad yy
                eval0.template values<0, false, false>(
                  hessians_quad + n_q_points, values_dofs + 2 * dofs_stride);
                // grad xy
                if ((evaluation_flag & EvaluationFlags::gradients))
                  eval0.template gradients<0, false, true>(
                    hessians_quad + 2 * n_q_points, values_dofs + dofs_stride);
                else
                  eval0.template gradients<0, false, false>(
                    hessians_quad + 2 * n_q_points, values_dofs + dofs_stride);
                break;
              default:
                AssertThrow(false, ExcNotImplemented());
            }
        }
    }
  };


  template <int dim, int fe_degree, typename Number>
  struct FEFaceNormalEvaluationImpl
  {
    template <bool do_evaluate, bool add_into_output>
    static void
    interpolate(const unsigned int                            n_components,
                const EvaluationFlags::EvaluationFlags        flags,
                const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                const Number *                                input,
                Number *                                      output,
                const unsigned int                            face_no)
    {
      Assert(static_cast<unsigned int>(fe_degree) ==
                 shape_info.data.front().fe_degree ||
               fe_degree == -1,
             ExcInternalError());
      if (shape_info.element_type == MatrixFreeFunctions::tensor_raviart_thomas)
        interpolate_generic_raviart_thomas<do_evaluate, add_into_output>(
          n_components, input, output, flags, face_no, shape_info);
      else
        interpolate_generic<do_evaluate, add_into_output>(
          n_components,
          input,
          output,
          flags,
          face_no,
          shape_info.data.front().fe_degree + 1,
          shape_info.data.front().shape_data_on_face,
          shape_info.dofs_per_component_on_cell,
          3 * shape_info.dofs_per_component_on_face);
    }

    /**
     * Interpolate the values on the cell quadrature points onto a face.
     */
    template <bool do_evaluate, bool add_into_output>
    static void
    interpolate_quadrature(
      const unsigned int                            n_components,
      const EvaluationFlags::EvaluationFlags        flags,
      const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
      const Number *                                input,
      Number *                                      output,
      const unsigned int                            face_no)
    {
      Assert(static_cast<unsigned int>(fe_degree + 1) ==
                 shape_info.data.front().n_q_points_1d ||
               fe_degree == -1,
             ExcInternalError());

      interpolate_generic<do_evaluate, add_into_output>(
        n_components,
        input,
        output,
        flags,
        face_no,
        shape_info.data.front().quadrature.size(),
        shape_info.data.front().quadrature_data_on_face,
        shape_info.n_q_points,
        shape_info.n_q_points_face);
    }

  private:
    template <bool do_evaluate, bool add_into_output, int face_direction = 0>
    static void
    interpolate_generic(const unsigned int                     n_components,
                        const Number *                         input,
                        Number *                               output,
                        const EvaluationFlags::EvaluationFlags flag,
                        const unsigned int                     face_no,
                        const unsigned int                     n_points_1d,
                        const std::array<AlignedVector<Number>, 2> &shape_data,
                        const unsigned int dofs_per_component_on_cell,
                        const unsigned int dofs_per_component_on_face)
    {
      if (face_direction == face_no / 2)
        {
          EvaluatorTensorProduct<evaluate_general,
                                 dim,
                                 fe_degree + 1,
                                 0,
                                 Number>
            evalf(shape_data[face_no % 2],
                  AlignedVector<Number>(),
                  AlignedVector<Number>(),
                  n_points_1d,
                  0);

          const unsigned int in_stride  = do_evaluate ?
                                            dofs_per_component_on_cell :
                                            dofs_per_component_on_face;
          const unsigned int out_stride = do_evaluate ?
                                            dofs_per_component_on_face :
                                            dofs_per_component_on_cell;

          for (unsigned int c = 0; c < n_components; ++c)
            {
              if (flag & EvaluationFlags::hessians)
                evalf.template apply_face<face_direction,
                                          do_evaluate,
                                          add_into_output,
                                          2>(input, output);
              else if (flag & EvaluationFlags::gradients)
                evalf.template apply_face<face_direction,
                                          do_evaluate,
                                          add_into_output,
                                          1>(input, output);
              else
                evalf.template apply_face<face_direction,
                                          do_evaluate,
                                          add_into_output,
                                          0>(input, output);
              input += in_stride;
              output += out_stride;
            }
        }
      else if (face_direction < dim)
        {
          interpolate_generic<do_evaluate,
                              add_into_output,
                              std::min(face_direction + 1, dim - 1)>(
            n_components,
            input,
            output,
            flag,
            face_no,
            n_points_1d,
            shape_data,
            dofs_per_component_on_cell,
            dofs_per_component_on_face);
        }
    }

    template <typename EvalType>
    static EvalType
    create_evaluator_tensor_product(
      const MatrixFreeFunctions::UnivariateShapeData<Number> &data,
      const unsigned int                                      face_no)
    {
      return EvalType(data.shape_data_on_face[face_no % 2],
                      AlignedVector<Number>(),
                      AlignedVector<Number>());
    }

    template <bool do_evaluate,
              bool add_into_output,
              int  face_direction = 0,
              int  max_derivative = 0>
    static void
    interpolate_generic_raviart_thomas(
      const unsigned int                            n_components,
      const Number *                                input,
      Number *                                      output,
      const EvaluationFlags::EvaluationFlags        flag,
      const unsigned int                            face_no,
      const MatrixFreeFunctions::ShapeInfo<Number> &shape_info)
    {
      if (dim == 1)
        {
          // This should never happen since the FE_RaviartThomasNodal is not
          // defined for dim = 1. It prevents compiler warnings of infinite
          // recursion.
          Assert(false, ExcInternalError());
          return;
        }

      bool increase_max_der = false;
      if ((flag & EvaluationFlags::hessians && max_derivative < 2) ||
          (flag & EvaluationFlags::gradients && max_derivative < 1))
        increase_max_der = true;

      if (face_direction == face_no / 2 && !increase_max_der)
        {
          interpolate_generic_raviart_thomas_apply_face<do_evaluate,
                                                        add_into_output,
                                                        face_direction,
                                                        max_derivative>(
            shape_info, face_no, input, output);
        }
      else if (face_direction == face_no / 2)
        {
          // Only increase max_derivative
          interpolate_generic_raviart_thomas<do_evaluate,
                                             add_into_output,
                                             face_direction,
                                             std::min(max_derivative + 1, 2)>(
            n_components, input, output, flag, face_no, shape_info);
        }
      else if (face_direction < dim)
        {
          if (increase_max_der)
            {
              interpolate_generic_raviart_thomas<
                do_evaluate,
                add_into_output,
                std::min(face_direction + 1, dim - 1),
                std::min(max_derivative + 1, 2)>(
                n_components, input, output, flag, face_no, shape_info);
            }
          else
            {
              interpolate_generic_raviart_thomas<do_evaluate,
                                                 add_into_output,
                                                 std::min(face_direction + 1,
                                                          dim - 1),
                                                 max_derivative>(
                n_components, input, output, flag, face_no, shape_info);
            }
        }
    }

    /* Help function for interpolate_generic_raviart_thomas */
    template <bool do_evaluate,
              bool add_into_output,
              int  face_direction,
              int  max_derivative>
    static inline void
    interpolate_generic_raviart_thomas_apply_face(
      const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
      const unsigned int                            face_no,
      const Number *                                input,
      Number *                                      output)
    {
      // These types are evaluators in either normal or tangential direction
      // depending on the face direction, with different normal directions for
      // the different components.
      using Evalf0 = typename std::conditional<
        face_direction == 0,
        EvaluatorTensorProductAnisotropic<evaluate_raviart_thomas,
                                          dim,
                                          (fe_degree == -1) ? 1 : fe_degree + 1,
                                          0,
                                          Number,
                                          0>,
        EvaluatorTensorProductAnisotropic<evaluate_raviart_thomas,
                                          dim,
                                          (fe_degree == -1) ? 1 : fe_degree,
                                          0,
                                          Number,
                                          0>>::type;
      using Evalf1 = typename std::conditional<
        face_direction == 1,
        EvaluatorTensorProductAnisotropic<evaluate_raviart_thomas,
                                          dim,
                                          (fe_degree == -1) ? 1 : fe_degree + 1,
                                          0,
                                          Number,
                                          1>,
        EvaluatorTensorProductAnisotropic<evaluate_raviart_thomas,
                                          dim,
                                          (fe_degree == -1) ? 1 : fe_degree,
                                          0,
                                          Number,
                                          1>>::type;
      using Evalf2 = typename std::conditional<
        face_direction == 2,
        EvaluatorTensorProductAnisotropic<evaluate_raviart_thomas,
                                          dim,
                                          (fe_degree == -1) ? 1 : fe_degree + 1,
                                          0,
                                          Number,
                                          2>,
        EvaluatorTensorProductAnisotropic<evaluate_raviart_thomas,
                                          dim,
                                          (fe_degree == -1) ? 1 : fe_degree,
                                          0,
                                          Number,
                                          2>>::type;

      Evalf0 evalf0 =
        create_evaluator_tensor_product<Evalf0>((face_direction == 0) ?
                                                  shape_info.data[0] :
                                                  shape_info.data[1],
                                                face_no);
      Evalf1 evalf1 =
        create_evaluator_tensor_product<Evalf1>((face_direction == 1) ?
                                                  shape_info.data[0] :
                                                  shape_info.data[1],
                                                face_no);
      Evalf2 evalf2 =
        create_evaluator_tensor_product<Evalf2>((face_direction == 2) ?
                                                  shape_info.data[0] :
                                                  shape_info.data[1],
                                                face_no);

      const unsigned int dofs_per_component_on_cell =
        shape_info.dofs_per_component_on_cell;
      const unsigned int dofs_per_component_on_face =
        3 * shape_info.dofs_per_component_on_face;

      // NOTE! dofs_per_component_on_face is in the tangent direction,
      // i.e (fe.degree+1)*fe.degree. Normal faces are only
      // fe.degree*fe.degree
      const unsigned int in_stride =
        do_evaluate ? dofs_per_component_on_cell : dofs_per_component_on_face;
      const unsigned int out_stride =
        do_evaluate ? dofs_per_component_on_face : dofs_per_component_on_cell;

      const unsigned int in_stride_after_normal =
        do_evaluate ?
          dofs_per_component_on_cell :
          dofs_per_component_on_face - 3 * Utilities::pow(fe_degree, dim - 2);
      const unsigned int out_stride_after_normal =
        do_evaluate ?
          dofs_per_component_on_face - 3 * Utilities::pow(fe_degree, dim - 2) :
          dofs_per_component_on_cell;

      evalf0.template apply_face<face_direction,
                                 do_evaluate,
                                 add_into_output,
                                 max_derivative>(input, output);
      // stride to next component
      input += (face_direction == 0) ? in_stride_after_normal : in_stride;
      output += (face_direction == 0) ? out_stride_after_normal : out_stride;

      evalf1.template apply_face<face_direction,
                                 do_evaluate,
                                 add_into_output,
                                 max_derivative>(input, output);

      if (dim == 3)
        {
          // stride to next component
          input += (face_direction == 1) ? in_stride_after_normal : in_stride;
          output +=
            (face_direction == 1) ? out_stride_after_normal : out_stride;

          evalf2.template apply_face<face_direction,
                                     do_evaluate,
                                     add_into_output,
                                     max_derivative>(input, output);
        }
    }
  };



  // internal helper function for reading data; base version of different types
  template <typename VectorizedArrayType, typename Number2>
  void
  do_vectorized_read(const Number2 *src_ptr, VectorizedArrayType &dst)
  {
    for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
      dst[v] = src_ptr[v];
  }



  // internal helper function for reading data; specialized version where we
  // can use a dedicated load function
  template <typename Number, std::size_t width>
  void
  do_vectorized_read(const Number *src_ptr, VectorizedArray<Number, width> &dst)
  {
    dst.load(src_ptr);
  }



  // internal helper function for reading data; base version of different types
  template <typename VectorizedArrayType, typename Number2>
  void
  do_vectorized_gather(const Number2 *      src_ptr,
                       const unsigned int * indices,
                       VectorizedArrayType &dst)
  {
    for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
      dst[v] = src_ptr[indices[v]];
  }



  // internal helper function for reading data; specialized version where we
  // can use a dedicated gather function
  template <typename Number, std::size_t width>
  void
  do_vectorized_gather(const Number *                  src_ptr,
                       const unsigned int *            indices,
                       VectorizedArray<Number, width> &dst)
  {
    dst.gather(src_ptr, indices);
  }



  // internal helper function for reading data; base version of different types
  template <typename VectorizedArrayType, typename Number2>
  void
  do_vectorized_add(const VectorizedArrayType src, Number2 *dst_ptr)
  {
    for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
      dst_ptr[v] += src[v];
  }



  // internal helper function for reading data; specialized version where we
  // can use a dedicated load function
  template <typename Number, std::size_t width>
  void
  do_vectorized_add(const VectorizedArray<Number, width> src, Number *dst_ptr)
  {
    VectorizedArray<Number, width> tmp;
    tmp.load(dst_ptr);
    (tmp + src).store(dst_ptr);
  }



  // internal helper function for reading data; base version of different types
  template <typename VectorizedArrayType, typename Number2>
  void
  do_vectorized_scatter_add(const VectorizedArrayType src,
                            const unsigned int *      indices,
                            Number2 *                 dst_ptr)
  {
    for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
      dst_ptr[indices[v]] += src[v];
  }



  // internal helper function for reading data; specialized version where we
  // can use a dedicated gather function
  template <typename Number, std::size_t width>
  void
  do_vectorized_scatter_add(const VectorizedArray<Number, width> src,
                            const unsigned int *                 indices,
                            Number *                             dst_ptr)
  {
#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS < 512
    for (unsigned int v = 0; v < width; ++v)
      dst_ptr[indices[v]] += src[v];
#else
    VectorizedArray<Number, width> tmp;
    tmp.gather(dst_ptr, indices);
    (tmp + src).scatter(indices, dst_ptr);
#endif
  }



  template <typename Number>
  void
  adjust_for_face_orientation(const unsigned int dim,
                              const unsigned int n_components,
                              const EvaluationFlags::EvaluationFlags flag,
                              const unsigned int *orientation,
                              const bool          integrate,
                              const std::size_t   n_q_points,
                              Number *            tmp_values,
                              Number *            values_quad,
                              Number *            gradients_quad,
                              Number *            hessians_quad)
  {
    for (unsigned int c = 0; c < n_components; ++c)
      {
        if (flag & EvaluationFlags::values)
          {
            if (integrate)
              for (unsigned int q = 0; q < n_q_points; ++q)
                tmp_values[q] = values_quad[c * n_q_points + orientation[q]];
            else
              for (unsigned int q = 0; q < n_q_points; ++q)
                tmp_values[orientation[q]] = values_quad[c * n_q_points + q];
            for (unsigned int q = 0; q < n_q_points; ++q)
              values_quad[c * n_q_points + q] = tmp_values[q];
          }
        if (flag & EvaluationFlags::gradients)
          for (unsigned int d = 0; d < dim; ++d)
            {
              if (integrate)
                for (unsigned int q = 0; q < n_q_points; ++q)
                  tmp_values[q] =
                    gradients_quad[(c * dim + d) * n_q_points + orientation[q]];
              else
                for (unsigned int q = 0; q < n_q_points; ++q)
                  tmp_values[orientation[q]] =
                    gradients_quad[(c * dim + d) * n_q_points + q];
              for (unsigned int q = 0; q < n_q_points; ++q)
                gradients_quad[(c * dim + d) * n_q_points + q] = tmp_values[q];
            }
        if (flag & EvaluationFlags::hessians)
          {
            const unsigned int hdim = (dim * (dim + 1)) / 2;
            for (unsigned int d = 0; d < hdim; ++d)
              {
                if (integrate)
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    tmp_values[q] = hessians_quad[(c * hdim + d) * n_q_points +
                                                  orientation[q]];
                else
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    tmp_values[orientation[q]] =
                      hessians_quad[(c * hdim + d) * n_q_points + q];
                for (unsigned int q = 0; q < n_q_points; ++q)
                  hessians_quad[(c * hdim + d) * n_q_points + q] =
                    tmp_values[q];
              }
          }
      }
  }



  template <typename Number, typename VectorizedArrayType>
  void
  adjust_for_face_orientation_per_lane(
    const unsigned int                     dim,
    const unsigned int                     n_components,
    const unsigned int                     v,
    const EvaluationFlags::EvaluationFlags flag,
    const unsigned int *                   orientation,
    const bool                             integrate,
    const std::size_t                      n_q_points,
    Number *                               tmp_values,
    VectorizedArrayType *                  values_quad,
    VectorizedArrayType *                  gradients_quad = nullptr,
    VectorizedArrayType *                  hessians_quad  = nullptr)
  {
    for (unsigned int c = 0; c < n_components; ++c)
      {
        if (flag & EvaluationFlags::values)
          {
            if (integrate)
              for (unsigned int q = 0; q < n_q_points; ++q)
                tmp_values[q] = values_quad[c * n_q_points + orientation[q]][v];
            else
              for (unsigned int q = 0; q < n_q_points; ++q)
                tmp_values[orientation[q]] = values_quad[c * n_q_points + q][v];
            for (unsigned int q = 0; q < n_q_points; ++q)
              values_quad[c * n_q_points + q][v] = tmp_values[q];
          }
        if (flag & EvaluationFlags::gradients)
          for (unsigned int d = 0; d < dim; ++d)
            {
              Assert(gradients_quad != nullptr, ExcInternalError());
              if (integrate)
                for (unsigned int q = 0; q < n_q_points; ++q)
                  tmp_values[q] = gradients_quad[(c * dim + d) * n_q_points +
                                                 orientation[q]][v];
              else
                for (unsigned int q = 0; q < n_q_points; ++q)
                  tmp_values[orientation[q]] =
                    gradients_quad[(c * dim + d) * n_q_points + q][v];
              for (unsigned int q = 0; q < n_q_points; ++q)
                gradients_quad[(c * dim + d) * n_q_points + q][v] =
                  tmp_values[q];
            }
        if (flag & EvaluationFlags::hessians)
          {
            Assert(hessians_quad != nullptr, ExcInternalError());
            const unsigned int hdim = (dim * (dim + 1)) / 2;
            for (unsigned int d = 0; d < hdim; ++d)
              {
                if (integrate)
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    tmp_values[q] = hessians_quad[(c * hdim + d) * n_q_points +
                                                  orientation[q]][v];
                else
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    tmp_values[orientation[q]] =
                      hessians_quad[(c * hdim + d) * n_q_points + q][v];
                for (unsigned int q = 0; q < n_q_points; ++q)
                  hessians_quad[(c * hdim + d) * n_q_points + q][v] =
                    tmp_values[q];
              }
          }
      }
  }



  template <int dim, typename Number>
  struct FEFaceEvaluationImplEvaluateSelector
  {
    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                     n_components,
        const EvaluationFlags::EvaluationFlags evaluation_flag,
        const Number *                         values_dofs,
        FEEvaluationData<dim, Number, true> &  fe_eval)
    {
      const auto &shape_info = fe_eval.get_shape_info();
      const auto &shape_data = shape_info.data.front();

      if (shape_info.element_type == MatrixFreeFunctions::tensor_none)
        {
          Assert((fe_eval.get_dof_access_index() ==
                    MatrixFreeFunctions::DoFInfo::dof_access_cell &&
                  fe_eval.is_interior_face() == false) == false,
                 ExcNotImplemented());

          const unsigned int face_no          = fe_eval.get_face_no();
          const unsigned int face_orientation = fe_eval.get_face_orientation();
          const std::size_t  n_dofs     = shape_info.dofs_per_component_on_cell;
          const std::size_t  n_q_points = shape_info.n_q_points_faces[face_no];

          using Eval =
            EvaluatorTensorProduct<evaluate_general, 1, 0, 0, Number, Number>;

          if (evaluation_flag & EvaluationFlags::values)
            {
              const auto shape_values =
                &shape_data.shape_values_face(face_no, face_orientation, 0);

              auto values_quad_ptr        = fe_eval.begin_values();
              auto values_dofs_actual_ptr = values_dofs;

              Eval eval(shape_values, nullptr, nullptr, n_dofs, n_q_points);
              for (unsigned int c = 0; c < n_components; ++c)
                {
                  eval.template values<0, true, false>(values_dofs_actual_ptr,
                                                       values_quad_ptr);

                  values_quad_ptr += n_q_points;
                  values_dofs_actual_ptr += n_dofs;
                }
            }

          if (evaluation_flag & EvaluationFlags::gradients)
            {
              auto gradients_quad_ptr     = fe_eval.begin_gradients();
              auto values_dofs_actual_ptr = values_dofs;

              std::array<const Number *, dim> shape_gradients;
              for (unsigned int d = 0; d < dim; ++d)
                shape_gradients[d] = &shape_data.shape_gradients_face(
                  face_no, face_orientation, d, 0);

              for (unsigned int c = 0; c < n_components; ++c)
                {
                  for (unsigned int d = 0; d < dim; ++d)
                    {
                      Eval eval(nullptr,
                                shape_gradients[d],
                                nullptr,
                                n_dofs,
                                n_q_points);

                      eval.template gradients<0, true, false>(
                        values_dofs_actual_ptr, gradients_quad_ptr);

                      gradients_quad_ptr += n_q_points;
                    }
                  values_dofs_actual_ptr += n_dofs;
                }
            }

          Assert(!(evaluation_flag & EvaluationFlags::hessians),
                 ExcNotImplemented());

          return true;
        }

      const unsigned int dofs_per_face =
        fe_degree > -1 ? Utilities::pow(fe_degree + 1, dim - 1) :
                         Utilities::pow(shape_data.fe_degree + 1, dim - 1);

      // Note: we always keep storage of values, 1st and 2nd derivatives in an
      // array, so reserve space for all three here
      Number *temp         = fe_eval.get_scratch_data().begin();
      Number *scratch_data = temp + 3 * n_components * dofs_per_face;

      bool use_vectorization = true;

      if (fe_eval.get_dof_access_index() ==
            MatrixFreeFunctions::DoFInfo::dof_access_cell &&
          fe_eval.is_interior_face() == false) // exterior faces in the ECL loop
        use_vectorization =
          fe_eval.get_cell_ids()[0] != numbers::invalid_unsigned_int &&
          std::all_of(fe_eval.get_cell_ids().begin() + 1,
                      fe_eval.get_cell_ids().end(),
                      [&](const auto &v) {
                        return v == fe_eval.get_cell_ids()[0] ||
                               v == numbers::invalid_unsigned_int;
                      });

      if (use_vectorization == false)
        {
          for (unsigned int v = 0; v < Number::size(); ++v)
            {
              // the loop breaks once an invalid_unsigned_int is hit for
              // all cases except the exterior faces in the ECL loop (where
              // some faces might be at the boundaries but others not)
              if (fe_eval.get_cell_ids()[v] == numbers::invalid_unsigned_int)
                {
                  for (unsigned int i = 0; i < 3 * n_components * dofs_per_face;
                       ++i)
                    temp[i][v] = 0;
                  continue;
                }

              FEFaceNormalEvaluationImpl<dim, fe_degree, Number>::
                template interpolate<true, false>(n_components,
                                                  evaluation_flag,
                                                  shape_info,
                                                  values_dofs,
                                                  scratch_data,
                                                  fe_eval.get_face_no(v));

              for (unsigned int i = 0; i < 3 * n_components * dofs_per_face;
                   ++i)
                temp[i][v] = scratch_data[i][v];
            }
        }
      else
        FEFaceNormalEvaluationImpl<dim, fe_degree, Number>::
          template interpolate<true, false>(n_components,
                                            evaluation_flag,
                                            shape_info,
                                            values_dofs,
                                            temp,
                                            fe_eval.get_face_no());

      const unsigned int     subface_index = fe_eval.get_subface_index();
      constexpr unsigned int n_q_points_1d_actual =
        fe_degree > -1 ? n_q_points_1d : 0;

      if (fe_degree >= 1 &&
          shape_info.element_type == MatrixFreeFunctions::tensor_raviart_thomas)
        {
          FEFaceEvaluationImplRaviartThomas<dim,
                                            (fe_degree == -1) ? 1 : fe_degree,
                                            (n_q_points_1d < 1) ? 1 :
                                                                  n_q_points_1d,
                                            Number>::
            template evaluate_or_integrate_in_face<false>(
              evaluation_flag,
              temp,
              fe_eval,
              scratch_data,
              subface_index,
              fe_eval.get_face_no());
        }
      else if (fe_degree > -1 &&
               subface_index >= GeometryInfo<dim>::max_children_per_cell &&
               shape_info.element_type <= MatrixFreeFunctions::tensor_symmetric)
        FEFaceEvaluationImpl<true,
                             dim,
                             fe_degree,
                             n_q_points_1d_actual,
                             Number>::evaluate_in_face(n_components,
                                                       evaluation_flag,
                                                       shape_data,
                                                       temp,
                                                       fe_eval.begin_values(),
                                                       fe_eval
                                                         .begin_gradients(),
                                                       fe_eval.begin_hessians(),
                                                       scratch_data,
                                                       subface_index);
      else
        FEFaceEvaluationImpl<false,
                             dim,
                             fe_degree,
                             n_q_points_1d_actual,
                             Number>::evaluate_in_face(n_components,
                                                       evaluation_flag,
                                                       shape_data,
                                                       temp,
                                                       fe_eval.begin_values(),
                                                       fe_eval
                                                         .begin_gradients(),
                                                       fe_eval.begin_hessians(),
                                                       scratch_data,
                                                       subface_index);

      if (use_vectorization == false)
        {
          for (unsigned int v = 0; v < Number::size(); ++v)
            {
              // the loop breaks once an invalid_unsigned_int is hit for
              // all cases except the exterior faces in the ECL loop (where
              // some faces might be at the boundaries but others not)
              if (fe_eval.get_cell_ids()[v] == numbers::invalid_unsigned_int)
                continue;

              if (fe_eval.get_face_orientation(v) != 0)
                adjust_for_face_orientation_per_lane(
                  dim,
                  n_components,
                  v,
                  evaluation_flag,
                  &fe_eval.get_shape_info().face_orientations_quad(
                    fe_eval.get_face_orientation(v), 0),
                  false,
                  shape_info.n_q_points_face,
                  &temp[0][0],
                  fe_eval.begin_values(),
                  fe_eval.begin_gradients(),
                  fe_eval.begin_hessians());
            }
        }
      else if (fe_eval.get_face_orientation() != 0)
        adjust_for_face_orientation(
          dim,
          n_components,
          evaluation_flag,
          &fe_eval.get_shape_info().face_orientations_quad(
            fe_eval.get_face_orientation(), 0),
          false,
          shape_info.n_q_points_face,
          temp,
          fe_eval.begin_values(),
          fe_eval.begin_gradients(),
          fe_eval.begin_hessians());

      return false;
    }
  };



  template <int dim, typename Number>
  struct FEFaceEvaluationImplIntegrateSelector
  {
    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                     n_components,
        const EvaluationFlags::EvaluationFlags integration_flag,
        Number *                               values_dofs,
        FEEvaluationData<dim, Number, true> &  fe_eval)
    {
      const auto &shape_info = fe_eval.get_shape_info();
      const auto &shape_data = shape_info.data.front();

      if (shape_info.element_type == MatrixFreeFunctions::tensor_none)
        {
          Assert((fe_eval.get_dof_access_index() ==
                    MatrixFreeFunctions::DoFInfo::dof_access_cell &&
                  fe_eval.is_interior_face() == false) == false,
                 ExcNotImplemented());

          const unsigned int face_no          = fe_eval.get_face_no();
          const unsigned int face_orientation = fe_eval.get_face_orientation();
          const std::size_t  n_dofs     = shape_info.dofs_per_component_on_cell;
          const std::size_t  n_q_points = shape_info.n_q_points_faces[face_no];

          using Eval =
            EvaluatorTensorProduct<evaluate_general, 1, 0, 0, Number, Number>;

          if (integration_flag & EvaluationFlags::values)
            {
              const auto shape_values =
                &shape_data.shape_values_face(face_no, face_orientation, 0);

              auto values_quad_ptr        = fe_eval.begin_values();
              auto values_dofs_actual_ptr = values_dofs;

              Eval eval(shape_values, nullptr, nullptr, n_dofs, n_q_points);
              for (unsigned int c = 0; c < n_components; ++c)
                {
                  eval.template values<0, false, false>(values_quad_ptr,
                                                        values_dofs_actual_ptr);

                  values_quad_ptr += n_q_points;
                  values_dofs_actual_ptr += n_dofs;
                }
            }

          if (integration_flag & EvaluationFlags::gradients)
            {
              auto gradients_quad_ptr     = fe_eval.begin_gradients();
              auto values_dofs_actual_ptr = values_dofs;

              std::array<const Number *, dim> shape_gradients;
              for (unsigned int d = 0; d < dim; ++d)
                shape_gradients[d] = &shape_data.shape_gradients_face(
                  face_no, face_orientation, d, 0);

              for (unsigned int c = 0; c < n_components; ++c)
                {
                  for (unsigned int d = 0; d < dim; ++d)
                    {
                      Eval eval(nullptr,
                                shape_gradients[d],
                                nullptr,
                                n_dofs,
                                n_q_points);

                      if (!(integration_flag & EvaluationFlags::values) &&
                          d == 0)
                        eval.template gradients<0, false, false>(
                          gradients_quad_ptr, values_dofs_actual_ptr);
                      else
                        eval.template gradients<0, false, true>(
                          gradients_quad_ptr, values_dofs_actual_ptr);

                      gradients_quad_ptr += n_q_points;
                    }
                  values_dofs_actual_ptr += n_dofs;
                }
            }

          Assert(!(integration_flag & EvaluationFlags::hessians),
                 ExcNotImplemented());

          return true;
        }

      const unsigned int dofs_per_face =
        fe_degree > -1 ? Utilities::pow(fe_degree + 1, dim - 1) :
                         Utilities::pow(shape_data.fe_degree + 1, dim - 1);

      Number *temp         = fe_eval.get_scratch_data().begin();
      Number *scratch_data = temp + 3 * n_components * dofs_per_face;

      bool use_vectorization = true;

      if (fe_eval.get_dof_access_index() ==
            MatrixFreeFunctions::DoFInfo::dof_access_cell &&
          fe_eval.is_interior_face() == false) // exterior faces in the ECL loop
        use_vectorization =
          fe_eval.get_cell_ids()[0] != numbers::invalid_unsigned_int &&
          std::all_of(fe_eval.get_cell_ids().begin() + 1,
                      fe_eval.get_cell_ids().end(),
                      [&](const auto &v) {
                        return v == fe_eval.get_cell_ids()[0] ||
                               v == numbers::invalid_unsigned_int;
                      });

      if (use_vectorization == false)
        {
          for (unsigned int v = 0; v < Number::size(); ++v)
            {
              // the loop breaks once an invalid_unsigned_int is hit for
              // all cases except the exterior faces in the ECL loop (where
              // some faces might be at the boundaries but others not)
              if (fe_eval.get_cell_ids()[v] == numbers::invalid_unsigned_int)
                continue;

              if (fe_eval.get_face_orientation(v) != 0)
                adjust_for_face_orientation_per_lane(
                  dim,
                  n_components,
                  v,
                  integration_flag,
                  &fe_eval.get_shape_info().face_orientations_quad(
                    fe_eval.get_face_orientation(v), 0),
                  true,
                  shape_info.n_q_points_face,
                  &temp[0][0],
                  fe_eval.begin_values(),
                  fe_eval.begin_gradients(),
                  fe_eval.begin_hessians());
            }
        }
      else if (fe_eval.get_face_orientation() != 0)
        adjust_for_face_orientation(
          dim,
          n_components,
          integration_flag,
          &fe_eval.get_shape_info().face_orientations_quad(
            fe_eval.get_face_orientation(), 0),
          true,
          shape_info.n_q_points_face,
          temp,
          fe_eval.begin_values(),
          fe_eval.begin_gradients(),
          fe_eval.begin_hessians());

      const unsigned int n_q_points_1d_actual =
        fe_degree > -1 ? n_q_points_1d : 0;
      const unsigned int subface_index = fe_eval.get_subface_index();

      if (fe_degree >= 1 &&
          shape_info.element_type == MatrixFreeFunctions::tensor_raviart_thomas)
        {
          FEFaceEvaluationImplRaviartThomas<dim,
                                            (fe_degree == -1) ? 1 : fe_degree,
                                            (n_q_points_1d < 1) ? 1 :
                                                                  n_q_points_1d,
                                            Number>::
            template evaluate_or_integrate_in_face<true>(integration_flag,
                                                         temp,
                                                         fe_eval,
                                                         scratch_data,
                                                         subface_index,
                                                         fe_eval.get_face_no());
        }
      else if (fe_degree > -1 &&
               fe_eval.get_subface_index() >=
                 GeometryInfo<dim - 1>::max_children_per_cell &&
               shape_info.element_type <= MatrixFreeFunctions::tensor_symmetric)
        FEFaceEvaluationImpl<
          true,
          dim,
          fe_degree,
          n_q_points_1d_actual,
          Number>::integrate_in_face(n_components,
                                     integration_flag,
                                     shape_data,
                                     temp,
                                     fe_eval.begin_values(),
                                     fe_eval.begin_gradients(),
                                     fe_eval.begin_hessians(),
                                     scratch_data,
                                     subface_index);
      else
        FEFaceEvaluationImpl<
          false,
          dim,
          fe_degree,
          n_q_points_1d_actual,
          Number>::integrate_in_face(n_components,
                                     integration_flag,
                                     shape_data,
                                     temp,
                                     fe_eval.begin_values(),
                                     fe_eval.begin_gradients(),
                                     fe_eval.begin_hessians(),
                                     scratch_data,
                                     subface_index);

      if (use_vectorization == false)
        {
          for (unsigned int v = 0; v < Number::size(); ++v)
            {
              // the loop breaks once an invalid_unsigned_int is hit for
              // all cases except the exterior faces in the ECL loop (where
              // some faces might be at the boundaries but others not)
              if (fe_eval.get_cell_ids()[v] == numbers::invalid_unsigned_int)
                continue;

              FEFaceNormalEvaluationImpl<dim, fe_degree, Number>::
                template interpolate<false, false>(n_components,
                                                   integration_flag,
                                                   shape_info,
                                                   values_dofs,
                                                   scratch_data,
                                                   fe_eval.get_face_no(v));

              for (unsigned int i = 0; i < 3 * n_components * dofs_per_face;
                   ++i)
                temp[i][v] = scratch_data[i][v];
            }
        }
      else
        FEFaceNormalEvaluationImpl<dim, fe_degree, Number>::
          template interpolate<false, false>(n_components,
                                             integration_flag,
                                             shape_info,
                                             temp,
                                             values_dofs,
                                             fe_eval.get_face_no());
      return false;
    }
  };



  template <int n_face_orientations,
            typename Processor,
            typename EvaluationData,
            const bool check_face_orientations = false>
  void
  fe_face_evaluation_process_and_io(
    Processor &                            proc,
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags evaluation_flag,
    typename Processor::Number2_ *         global_vector_ptr,
    const std::vector<ArrayView<const typename Processor::Number2_>> *sm_ptr,
    const EvaluationData &                                            fe_eval,
    typename Processor::VectorizedArrayType_ *                        temp1)
  {
    constexpr int dim         = Processor::dim_;
    constexpr int fe_degree   = Processor::fe_degree_;
    using VectorizedArrayType = typename Processor::VectorizedArrayType_;
    constexpr int n_lanes     = VectorizedArrayType::size();

    using Number   = typename Processor::Number_;
    using Number2_ = typename Processor::Number2_;

    const auto &       shape_data = fe_eval.get_shape_info().data.front();
    constexpr bool     integrate  = Processor::do_integrate;
    const unsigned int face_no    = fe_eval.get_face_no();
    const auto &       dof_info   = fe_eval.get_dof_info();
    const unsigned int cell       = fe_eval.get_cell_or_face_batch_id();
    const MatrixFreeFunctions::DoFInfo::DoFAccessIndex dof_access_index =
      fe_eval.get_dof_access_index();
    AssertIndexRange(cell,
                     dof_info.index_storage_variants[dof_access_index].size());
    constexpr unsigned int dofs_per_face =
      Utilities::pow(fe_degree + 1, dim - 1);
    const unsigned int subface_index = fe_eval.get_subface_index();

    const unsigned int n_filled_lanes =
      dof_info.n_vectorization_lanes_filled[dof_access_index][cell];

    bool all_faces_are_same = n_filled_lanes == n_lanes;
    if (n_face_orientations == n_lanes)
      for (unsigned int v = 1; v < n_lanes; ++v)
        if (fe_eval.get_face_no(v) != fe_eval.get_face_no(0) ||
            fe_eval.get_face_orientation(v) != fe_eval.get_face_orientation(0))
          {
            all_faces_are_same = false;
            break;
          }

    // check for re-orientation ...
    std::array<const unsigned int *, n_face_orientations> orientation = {};

    if (dim == 3 && n_face_orientations == n_lanes && !all_faces_are_same &&
        fe_eval.is_interior_face() == 0)
      for (unsigned int v = 0; v < n_lanes; ++v)
        {
          // the loop breaks once an invalid_unsigned_int is hit for
          // all cases except the exterior faces in the ECL loop (where
          // some faces might be at the boundaries but others not)
          if (fe_eval.get_cell_ids()[v] == numbers::invalid_unsigned_int)
            continue;

          if (shape_data.nodal_at_cell_boundaries &&
              fe_eval.get_face_orientation(v) != 0)
            {
              // ... and in case we detect a re-orientation, go to the other
              // version of this function that actually allows for this
              if (subface_index == GeometryInfo<dim>::max_children_per_cell &&
                  check_face_orientations == false)
                {
                  fe_face_evaluation_process_and_io<n_face_orientations,
                                                    Processor,
                                                    EvaluationData,
                                                    true>(proc,
                                                          n_components,
                                                          evaluation_flag,
                                                          global_vector_ptr,
                                                          sm_ptr,
                                                          fe_eval,
                                                          temp1);
                  return;
                }
              orientation[v] = &fe_eval.get_shape_info().face_orientations_dofs(
                fe_eval.get_face_orientation(v), 0);
            }
        }
    else if (dim == 3 && fe_eval.get_face_orientation() != 0)
      {
        // go to the other version of this function
        if (subface_index == GeometryInfo<dim>::max_children_per_cell &&
            check_face_orientations == false)
          {
            fe_face_evaluation_process_and_io<n_face_orientations,
                                              Processor,
                                              EvaluationData,
                                              true>(proc,
                                                    n_components,
                                                    evaluation_flag,
                                                    global_vector_ptr,
                                                    sm_ptr,
                                                    fe_eval,
                                                    temp1);
            return;
          }
        for (unsigned int v = 0; v < n_face_orientations; ++v)
          orientation[v] = &fe_eval.get_shape_info().face_orientations_dofs(
            fe_eval.get_face_orientation(), 0);
      }

    // we know that the gradient weights for the Hermite case on the
    // right (side==1) are the negative from the value at the left
    // (side==0), so we only read out one of them.
    VectorizedArrayType grad_weight =
      shape_data
        .shape_data_on_face[0][fe_degree + (integrate ? (2 - face_no % 2) :
                                                        (1 + face_no % 2))];

    // face_to_cell_index_hermite
    std::array<const unsigned int *, n_face_orientations> index_array_hermite =
      {};
    if (fe_degree > 1 && (evaluation_flag & EvaluationFlags::gradients))
      {
        if (n_face_orientations == 1)
          index_array_hermite[0] =
            &fe_eval.get_shape_info().face_to_cell_index_hermite(face_no, 0);
        else
          {
            for (unsigned int v = 0; v < n_lanes; ++v)
              {
                if (fe_eval.get_cell_ids()[v] == numbers::invalid_unsigned_int)
                  continue;

                const auto face_no = fe_eval.get_face_no(v);

                grad_weight[v] =
                  shape_data.shape_data_on_face[0][fe_degree +
                                                   (integrate ?
                                                      (2 - (face_no % 2)) :
                                                      (1 + (face_no % 2)))][0];

                index_array_hermite[v] =
                  &fe_eval.get_shape_info().face_to_cell_index_hermite(face_no,
                                                                       0);
              }
          }
      }

    // face_to_cell_index_nodal
    std::array<const unsigned int *, n_face_orientations> index_array_nodal =
      {};
    if (shape_data.nodal_at_cell_boundaries == true)
      {
        if (n_face_orientations == 1)
          index_array_nodal[0] =
            &fe_eval.get_shape_info().face_to_cell_index_nodal(face_no, 0);
        else
          {
            for (unsigned int v = 0; v < n_lanes; ++v)
              {
                if (fe_eval.get_cell_ids()[v] == numbers::invalid_unsigned_int)
                  continue;

                const auto face_no = fe_eval.get_face_no(v);

                index_array_nodal[v] =
                  &fe_eval.get_shape_info().face_to_cell_index_nodal(face_no,
                                                                     0);
              }
          }
      }


    const auto reorientate = [&](const unsigned int v, const unsigned int i) {
      return (!check_face_orientations || orientation[v] == nullptr) ?
               i :
               orientation[v][i];
    };

    const unsigned int cell_index =
      dof_access_index == MatrixFreeFunctions::DoFInfo::dof_access_cell ?
        fe_eval.get_cell_ids()[0] :
        cell * n_lanes;
    const unsigned int *dof_indices =
      &dof_info.dof_indices_contiguous[dof_access_index][cell_index];

    for (unsigned int comp = 0; comp < n_components; ++comp)
      {
        const std::size_t index_offset =
          dof_info.component_dof_indices_offset
            [fe_eval.get_active_fe_index()]
            [fe_eval.get_first_selected_component()] +
          comp * Utilities::pow(fe_degree + 1, dim);

        // case 1: contiguous and interleaved indices
        if (n_face_orientations == 1 &&
            dof_info.index_storage_variants[dof_access_index][cell] ==
              MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
                interleaved_contiguous)
          {
            AssertDimension(
              dof_info.n_vectorization_lanes_filled[dof_access_index][cell],
              n_lanes);
            Number2_ *vector_ptr =
              global_vector_ptr + dof_indices[0] + index_offset * n_lanes;

            if (fe_degree > 1 && (evaluation_flag & EvaluationFlags::gradients))
              {
                for (unsigned int i = 0; i < dofs_per_face; ++i)
                  {
                    Assert(n_face_orientations == 1, ExcNotImplemented());

                    const unsigned int ind1 = index_array_hermite[0][2 * i];
                    const unsigned int ind2 = index_array_hermite[0][2 * i + 1];
                    const unsigned int i_   = reorientate(0, i);
                    proc.hermite_grad_vectorized(temp1[i_],
                                                 temp1[i_ + dofs_per_face],
                                                 vector_ptr + ind1 * n_lanes,
                                                 vector_ptr + ind2 * n_lanes,
                                                 grad_weight);
                  }
              }
            else
              {
                for (unsigned int i = 0; i < dofs_per_face; ++i)
                  {
                    Assert(n_face_orientations == 1, ExcNotImplemented());

                    const unsigned int i_  = reorientate(0, i);
                    const unsigned int ind = index_array_nodal[0][i];
                    proc.value_vectorized(temp1[i_],
                                          vector_ptr + ind * n_lanes);
                  }
              }
          }

        // case 2: contiguous and interleaved indices with fixed stride
        else if (n_face_orientations == 1 &&
                 dof_info.index_storage_variants[dof_access_index][cell] ==
                   MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
                     interleaved_contiguous_strided)
          {
            AssertDimension(
              dof_info.n_vectorization_lanes_filled[dof_access_index][cell],
              n_lanes);
            Number2_ *vector_ptr = global_vector_ptr + index_offset * n_lanes;
            if (fe_degree > 1 && (evaluation_flag & EvaluationFlags::gradients))
              {
                for (unsigned int i = 0; i < dofs_per_face; ++i)
                  {
                    Assert(n_face_orientations == 1, ExcNotImplemented());

                    const unsigned int i_ = reorientate(0, i);
                    const unsigned int ind1 =
                      index_array_hermite[0][2 * i] * n_lanes;
                    const unsigned int ind2 =
                      index_array_hermite[0][2 * i + 1] * n_lanes;
                    proc.hermite_grad_vectorized_indexed(
                      temp1[i_],
                      temp1[i_ + dofs_per_face],
                      vector_ptr + ind1,
                      vector_ptr + ind2,
                      grad_weight,
                      dof_indices,
                      dof_indices);
                  }
              }
            else
              {
                for (unsigned int i = 0; i < dofs_per_face; ++i)
                  {
                    Assert(n_face_orientations == 1, ExcNotImplemented());

                    const unsigned int i_  = reorientate(0, i);
                    const unsigned int ind = index_array_nodal[0][i] * n_lanes;
                    proc.value_vectorized_indexed(temp1[i_],
                                                  vector_ptr + ind,
                                                  dof_indices);
                  }
              }
          }

        // case 3: contiguous and interleaved indices with mixed stride
        else if (n_face_orientations == 1 &&
                 dof_info.index_storage_variants[dof_access_index][cell] ==
                   MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
                     interleaved_contiguous_mixed_strides)
          {
            const unsigned int *strides =
              &dof_info.dof_indices_interleave_strides[dof_access_index]
                                                      [cell * n_lanes];
            unsigned int indices[n_lanes];
            for (unsigned int v = 0; v < n_lanes; ++v)
              indices[v] = dof_indices[v] + index_offset * strides[v];
            const unsigned int n_filled_lanes =
              dof_info.n_vectorization_lanes_filled[dof_access_index][cell];

            if (fe_degree > 1 && (evaluation_flag & EvaluationFlags::gradients))
              {
                if (n_filled_lanes == n_lanes)
                  for (unsigned int i = 0; i < dofs_per_face; ++i)
                    {
                      Assert(n_face_orientations == 1, ExcNotImplemented());

                      const unsigned int i_ = reorientate(0, i);
                      unsigned int       ind1[n_lanes];
                      DEAL_II_OPENMP_SIMD_PRAGMA
                      for (unsigned int v = 0; v < n_lanes; ++v)
                        ind1[v] = indices[v] +
                                  index_array_hermite[0][2 * i] * strides[v];
                      unsigned int ind2[n_lanes];
                      DEAL_II_OPENMP_SIMD_PRAGMA
                      for (unsigned int v = 0; v < n_lanes; ++v)
                        ind2[v] =
                          indices[v] +
                          // TODO
                          index_array_hermite[0][2 * i + 1] * strides[v];
                      proc.hermite_grad_vectorized_indexed(
                        temp1[i_],
                        temp1[i_ + dofs_per_face],
                        global_vector_ptr,
                        global_vector_ptr,
                        grad_weight,
                        ind1,
                        ind2);
                    }
                else
                  {
                    if (integrate == false)
                      for (unsigned int i = 0; i < 2 * dofs_per_face; ++i)
                        temp1[i] = VectorizedArrayType();

                    for (unsigned int v = 0; v < n_filled_lanes; ++v)
                      for (unsigned int i = 0; i < dofs_per_face; ++i)
                        {
                          const unsigned int i_ =
                            reorientate(n_face_orientations == 1 ? 0 : v, i);
                          proc.hermite_grad(
                            temp1[i_][v],
                            temp1[i_ + dofs_per_face][v],
                            global_vector_ptr
                              [indices[v] +
                               index_array_hermite
                                   [n_face_orientations == 1 ? 0 : v][2 * i] *
                                 strides[v]],
                            global_vector_ptr
                              [indices[v] +
                               index_array_hermite[n_face_orientations == 1 ?
                                                     0 :
                                                     v][2 * i + 1] *
                                 strides[v]],
                            grad_weight[n_face_orientations == 1 ? 0 : v]);
                        }
                  }
              }
            else
              {
                if (n_filled_lanes == n_lanes)
                  for (unsigned int i = 0; i < dofs_per_face; ++i)
                    {
                      Assert(n_face_orientations == 1, ExcInternalError());
                      unsigned int ind[n_lanes];
                      DEAL_II_OPENMP_SIMD_PRAGMA
                      for (unsigned int v = 0; v < n_lanes; ++v)
                        ind[v] =
                          indices[v] + index_array_nodal[0][i] * strides[v];
                      const unsigned int i_ = reorientate(0, i);
                      proc.value_vectorized_indexed(temp1[i_],
                                                    global_vector_ptr,
                                                    ind);
                    }
                else
                  {
                    if (integrate == false)
                      for (unsigned int i = 0; i < dofs_per_face; ++i)
                        temp1[i] = VectorizedArrayType();

                    for (unsigned int v = 0; v < n_filled_lanes; ++v)
                      for (unsigned int i = 0; i < dofs_per_face; ++i)
                        proc.value(
                          temp1[reorientate(n_face_orientations == 1 ? 0 : v,
                                            i)][v],
                          global_vector_ptr
                            [indices[v] +
                             index_array_nodal[n_face_orientations == 1 ? 0 : v]
                                              [i] *
                               strides[v]]);
                  }
              }
          }

        // case 4: contiguous indices without interleaving
        else if (n_face_orientations > 1 ||
                 dof_info.index_storage_variants[dof_access_index][cell] ==
                   MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
                     contiguous)
          {
            Number2_ *vector_ptr = global_vector_ptr + index_offset;

            const bool vectorization_possible =
              all_faces_are_same && (sm_ptr == nullptr);

            std::array<Number2_ *, n_lanes>   vector_ptrs;
            std::array<unsigned int, n_lanes> reordered_indices;

            if (vectorization_possible == false)
              {
                vector_ptrs = {};
                if (n_face_orientations == 1)
                  {
                    for (unsigned int v = 0; v < n_filled_lanes; ++v)
                      if (sm_ptr == nullptr)
                        {
                          vector_ptrs[v] = vector_ptr + dof_indices[v];
                        }
                      else
                        {
                          const auto &temp =
                            dof_info
                              .dof_indices_contiguous_sm[dof_access_index]
                                                        [cell * n_lanes + v];
                          vector_ptrs[v] = const_cast<Number2_ *>(
                            sm_ptr->operator[](temp.first).data() +
                            temp.second + index_offset);
                        }
                  }
                else if (n_face_orientations == n_lanes)
                  {
                    const auto &cells = fe_eval.get_cell_ids();
                    for (unsigned int v = 0; v < n_lanes; ++v)
                      if (cells[v] != numbers::invalid_unsigned_int)
                        {
                          if (sm_ptr == nullptr)
                            {
                              vector_ptrs[v] =
                                vector_ptr +
                                dof_info
                                  .dof_indices_contiguous[dof_access_index]
                                                         [cells[v]];
                            }
                          else
                            {
                              const auto &temp =
                                dof_info
                                  .dof_indices_contiguous_sm[dof_access_index]
                                                            [cells[v]];
                              vector_ptrs[v] = const_cast<Number2_ *>(
                                sm_ptr->operator[](temp.first).data() +
                                temp.second + index_offset);
                            }
                        }
                  }
                else
                  {
                    Assert(false, ExcNotImplemented());
                  }
              }
            else if (n_face_orientations == n_lanes)
              {
                for (unsigned int v = 0; v < n_lanes; ++v)
                  reordered_indices[v] =
                    dof_info.dof_indices_contiguous[dof_access_index]
                                                   [fe_eval.get_cell_ids()[v]];
                dof_indices = reordered_indices.data();
              }

            if (fe_degree > 1 && (evaluation_flag & EvaluationFlags::gradients))
              {
                if (vectorization_possible)
                  for (unsigned int i = 0; i < dofs_per_face; ++i)
                    {
                      const unsigned int ind1 = index_array_hermite[0][2 * i];
                      const unsigned int ind2 =
                        index_array_hermite[0][2 * i + 1];
                      const unsigned int i_ = reorientate(0, i);

                      proc.hermite_grad_vectorized_indexed(
                        temp1[i_],
                        temp1[i_ + dofs_per_face],
                        vector_ptr + ind1,
                        vector_ptr + ind2,
                        grad_weight,
                        dof_indices,
                        dof_indices);
                    }
                else if (n_face_orientations == 1)
                  for (unsigned int i = 0; i < dofs_per_face; ++i)
                    {
                      const unsigned int ind1 = index_array_hermite[0][2 * i];
                      const unsigned int ind2 =
                        index_array_hermite[0][2 * i + 1];
                      const unsigned int i_ = reorientate(0, i);

                      for (unsigned int v = 0; v < n_filled_lanes; ++v)
                        proc.hermite_grad(temp1[i_][v],
                                          temp1[i_ + dofs_per_face][v],
                                          vector_ptrs[v][ind1],
                                          vector_ptrs[v][ind2],
                                          grad_weight[v]);

                      if (integrate == false)
                        for (unsigned int v = n_filled_lanes; v < n_lanes; ++v)
                          {
                            temp1[i][v]                 = 0.0;
                            temp1[i + dofs_per_face][v] = 0.0;
                          }
                    }
                else
                  {
                    if (integrate == false && n_filled_lanes < n_lanes)
                      for (unsigned int i = 0; i < dofs_per_face; ++i)
                        temp1[i] = temp1[i + dofs_per_face] = Number();

                    for (unsigned int v = 0; v < n_filled_lanes; ++v)
                      for (unsigned int i = 0; i < dofs_per_face; ++i)
                        proc.hermite_grad(
                          temp1[reorientate(v, i)][v],
                          temp1[reorientate(v, i) + dofs_per_face][v],
                          vector_ptrs[v][index_array_hermite[v][2 * i]],
                          vector_ptrs[v][index_array_hermite[v][2 * i + 1]],
                          grad_weight[v]);
                  }
              }
            else
              {
                if (vectorization_possible)
                  for (unsigned int i = 0; i < dofs_per_face; ++i)
                    {
                      const unsigned int ind = index_array_nodal[0][i];
                      const unsigned int i_  = reorientate(0, i);

                      proc.value_vectorized_indexed(temp1[i_],
                                                    vector_ptr + ind,
                                                    dof_indices);
                    }
                else if (n_face_orientations == 1)
                  for (unsigned int i = 0; i < dofs_per_face; ++i)
                    {
                      const unsigned int ind = index_array_nodal[0][i];
                      const unsigned int i_  = reorientate(0, i);

                      for (unsigned int v = 0; v < n_filled_lanes; ++v)
                        proc.value(temp1[i_][v], vector_ptrs[v][ind]);

                      if (integrate == false)
                        for (unsigned int v = n_filled_lanes; v < n_lanes; ++v)
                          temp1[i_][v] = 0.0;
                    }
                else
                  {
                    if (integrate == false && n_filled_lanes < n_lanes)
                      for (unsigned int i = 0; i < dofs_per_face; ++i)
                        temp1[i] = Number();

                    for (unsigned int v = 0; v < n_filled_lanes; ++v)
                      for (unsigned int i = 0; i < dofs_per_face; ++i)
                        proc.value(temp1[reorientate(v, i)][v],
                                   vector_ptrs[v][index_array_nodal[v][i]]);
                  }
              }
          }
        else
          {
            // We should not end up here, this should be caught by
            // FEFaceEvaluationImplGatherEvaluateSelector::supports()
            Assert(false, ExcInternalError());
          }
        temp1 += 3 * dofs_per_face;
      }
  }



  template <int dim, typename Number2, typename VectorizedArrayType>
  struct FEFaceEvaluationImplGatherEvaluateSelector
  {
    using Number = typename VectorizedArrayType::value_type;

    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                                n_components,
        const EvaluationFlags::EvaluationFlags            evaluation_flag,
        const Number2 *                                   src_ptr,
        const std::vector<ArrayView<const Number2>> *     sm_ptr,
        FEEvaluationData<dim, VectorizedArrayType, true> &fe_eval)
    {
      Assert(fe_degree > -1, ExcInternalError());
      Assert(fe_eval.get_shape_info().element_type <=
               MatrixFreeFunctions::tensor_symmetric,
             ExcInternalError());

      const unsigned int dofs_per_face = Utilities::pow(fe_degree + 1, dim - 1);

      VectorizedArrayType *temp = fe_eval.get_scratch_data().begin();
      VectorizedArrayType *scratch_data =
        temp + 3 * n_components * dofs_per_face;

      Processor<fe_degree> p;

      if (fe_eval.get_dof_access_index() ==
            MatrixFreeFunctions::DoFInfo::dof_access_cell &&
          fe_eval.is_interior_face() == false)
        fe_face_evaluation_process_and_io<VectorizedArrayType::size()>(
          p, n_components, evaluation_flag, src_ptr, sm_ptr, fe_eval, temp);
      else
        fe_face_evaluation_process_and_io<1>(
          p, n_components, evaluation_flag, src_ptr, sm_ptr, fe_eval, temp);

      const unsigned int subface_index = fe_eval.get_subface_index();

      if (subface_index >= GeometryInfo<dim>::max_children_per_cell)
        FEFaceEvaluationImpl<true,
                             dim,
                             fe_degree,
                             n_q_points_1d,
                             VectorizedArrayType>::
          evaluate_in_face(n_components,
                           evaluation_flag,
                           fe_eval.get_shape_info().data.front(),
                           temp,
                           fe_eval.begin_values(),
                           fe_eval.begin_gradients(),
                           fe_eval.begin_hessians(),
                           scratch_data,
                           subface_index);
      else
        FEFaceEvaluationImpl<false,
                             dim,
                             fe_degree,
                             n_q_points_1d,
                             VectorizedArrayType>::
          evaluate_in_face(n_components,
                           evaluation_flag,
                           fe_eval.get_shape_info().data.front(),
                           temp,
                           fe_eval.begin_values(),
                           fe_eval.begin_gradients(),
                           fe_eval.begin_hessians(),
                           scratch_data,
                           subface_index);

      // re-orientation for cases not possible with above algorithm
      if (subface_index < GeometryInfo<dim>::max_children_per_cell)
        {
          if (fe_eval.get_dof_access_index() ==
                MatrixFreeFunctions::DoFInfo::dof_access_cell &&
              fe_eval.is_interior_face() == false)
            {
              for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
                {
                  // the loop breaks once an invalid_unsigned_int is hit for
                  // all cases except the exterior faces in the ECL loop (where
                  // some faces might be at the boundaries but others not)
                  if (fe_eval.get_cell_ids()[v] ==
                      numbers::invalid_unsigned_int)
                    continue;

                  if (fe_eval.get_face_orientation(v) != 0)
                    adjust_for_face_orientation_per_lane(
                      dim,
                      n_components,
                      v,
                      evaluation_flag,
                      &fe_eval.get_shape_info().face_orientations_quad(
                        fe_eval.get_face_orientation(v), 0),
                      false,
                      Utilities::pow(n_q_points_1d, dim - 1),
                      &temp[0][0],
                      fe_eval.begin_values(),
                      fe_eval.begin_gradients(),
                      fe_eval.begin_hessians());
                }
            }
          else if (fe_eval.get_face_orientation() != 0)
            adjust_for_face_orientation(
              dim,
              n_components,
              evaluation_flag,
              &fe_eval.get_shape_info().face_orientations_quad(
                fe_eval.get_face_orientation(), 0),
              false,
              Utilities::pow(n_q_points_1d, dim - 1),
              temp,
              fe_eval.begin_values(),
              fe_eval.begin_gradients(),
              fe_eval.begin_hessians());
        }

      return false;
    }

    static bool
    supports(
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &shape_info,
      const Number2 *                                            vector_ptr,
      MatrixFreeFunctions::DoFInfo::IndexStorageVariants         storage)
    {
      const unsigned int fe_degree = shape_info.data.front().fe_degree;
      if (fe_degree < 1 || !shape_info.data.front().nodal_at_cell_boundaries ||
          (evaluation_flag & EvaluationFlags::gradients &&
           (fe_degree < 2 ||
            shape_info.data.front().element_type !=
              MatrixFreeFunctions::tensor_symmetric_hermite)) ||
          (evaluation_flag & EvaluationFlags::hessians) ||
          vector_ptr == nullptr ||
          shape_info.data.front().element_type >
            MatrixFreeFunctions::tensor_symmetric ||
          storage <
            MatrixFreeFunctions::DoFInfo::IndexStorageVariants::contiguous)
        return false;
      else
        return true;
    }

  private:
    template <int fe_degree>
    struct Processor
    {
      static const bool do_integrate = false;
      static const int  dim_         = dim;
      static const int  fe_degree_   = fe_degree;
      using VectorizedArrayType_     = VectorizedArrayType;
      using Number_                  = Number;
      using Number2_                 = const Number2;

      template <typename T0, typename T1, typename T2>
      void
      hermite_grad_vectorized(T0 &      temp_1,
                              T0 &      temp_2,
                              const T1  src_ptr_1,
                              const T1  src_ptr_2,
                              const T2 &grad_weight)
      {
        do_vectorized_read(src_ptr_1, temp_1);
        do_vectorized_read(src_ptr_2, temp_2);
        temp_2 = grad_weight * (temp_1 - temp_2);
      }

      template <typename T1, typename T2>
      void
      value_vectorized(T1 &temp, const T2 src_ptr)
      {
        do_vectorized_read(src_ptr, temp);
      }

      template <typename T0, typename T1, typename T2, typename T3>
      void
      hermite_grad_vectorized_indexed(T0 &      temp_1,
                                      T0 &      temp_2,
                                      const T1  src_ptr_1,
                                      const T1  src_ptr_2,
                                      const T2 &grad_weight,
                                      const T3 &indices_1,
                                      const T3 &indices_2)
      {
        do_vectorized_gather(src_ptr_1, indices_1, temp_1);
        do_vectorized_gather(src_ptr_2, indices_2, temp_2);
        temp_2 = grad_weight * (temp_1 - temp_2);
      }

      template <typename T0, typename T1, typename T2>
      void
      value_vectorized_indexed(T0 &temp, const T1 src_ptr, const T2 &indices)
      {
        do_vectorized_gather(src_ptr, indices, temp);
      }

      template <typename T0, typename T1, typename T2>
      void
      hermite_grad(T0 &      temp_1,
                   T0 &      temp_2,
                   const T1 &src_ptr_1,
                   const T1 &src_ptr_2,
                   const T2 &grad_weight)
      {
        // case 3a)
        temp_1 = src_ptr_1;
        temp_2 = grad_weight * (temp_1 - src_ptr_2);
      }

      template <typename T1, typename T2>
      void
      value(T1 &temp, const T2 &src_ptr)
      {
        // case 3b)
        temp = src_ptr;
      }
    };
  };



  template <int dim, typename Number2, typename VectorizedArrayType>
  struct FEFaceEvaluationImplIntegrateScatterSelector
  {
    using Number = typename VectorizedArrayType::value_type;

    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                                n_components,
        const EvaluationFlags::EvaluationFlags            integration_flag,
        Number2 *                                         dst_ptr,
        const std::vector<ArrayView<const Number2>> *     sm_ptr,
        FEEvaluationData<dim, VectorizedArrayType, true> &fe_eval)
    {
      Assert(fe_degree > -1, ExcInternalError());
      Assert(fe_eval.get_shape_info().element_type <=
               MatrixFreeFunctions::tensor_symmetric,
             ExcInternalError());

      const unsigned int dofs_per_face = Utilities::pow(fe_degree + 1, dim - 1);

      VectorizedArrayType *temp = fe_eval.get_scratch_data().begin();
      VectorizedArrayType *scratch_data =
        temp + 3 * n_components * dofs_per_face;

      const unsigned int subface_index = fe_eval.get_subface_index();

      // re-orientation for cases not possible with the io function below
      if (subface_index < GeometryInfo<dim>::max_children_per_cell)
        {
          if (fe_eval.get_dof_access_index() ==
                MatrixFreeFunctions::DoFInfo::dof_access_cell &&
              fe_eval.is_interior_face() == false)
            for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
              {
                // the loop breaks once an invalid_unsigned_int is hit for
                // all cases except the exterior faces in the ECL loop (where
                // some faces might be at the boundaries but others not)
                if (fe_eval.get_cell_ids()[v] == numbers::invalid_unsigned_int)
                  continue;

                if (fe_eval.get_face_orientation(v) != 0)
                  adjust_for_face_orientation_per_lane(
                    dim,
                    n_components,
                    v,
                    integration_flag,
                    &fe_eval.get_shape_info().face_orientations_quad(
                      fe_eval.get_face_orientation(v), 0),
                    true,
                    Utilities::pow(n_q_points_1d, dim - 1),
                    &temp[0][0],
                    fe_eval.begin_values(),
                    fe_eval.begin_gradients(),
                    fe_eval.begin_hessians());
              }
          else if (fe_eval.get_face_orientation() != 0)
            adjust_for_face_orientation(
              dim,
              n_components,
              integration_flag,
              &fe_eval.get_shape_info().face_orientations_quad(
                fe_eval.get_face_orientation(), 0),
              true,
              Utilities::pow(n_q_points_1d, dim - 1),
              temp,
              fe_eval.begin_values(),
              fe_eval.begin_gradients(),
              fe_eval.begin_hessians());
        }

      if (fe_degree > -1 && fe_eval.get_subface_index() >=
                              GeometryInfo<dim - 1>::max_children_per_cell)
        FEFaceEvaluationImpl<true,
                             dim,
                             fe_degree,
                             n_q_points_1d,
                             VectorizedArrayType>::
          integrate_in_face(n_components,
                            integration_flag,
                            fe_eval.get_shape_info().data.front(),
                            temp,
                            fe_eval.begin_values(),
                            fe_eval.begin_gradients(),
                            fe_eval.begin_hessians(),
                            scratch_data,
                            subface_index);
      else
        FEFaceEvaluationImpl<false,
                             dim,
                             fe_degree,
                             n_q_points_1d,
                             VectorizedArrayType>::
          integrate_in_face(n_components,
                            integration_flag,
                            fe_eval.get_shape_info().data.front(),
                            temp,
                            fe_eval.begin_values(),
                            fe_eval.begin_gradients(),
                            fe_eval.begin_hessians(),
                            scratch_data,
                            subface_index);

      Processor<fe_degree> p;

      if (fe_eval.get_dof_access_index() ==
            MatrixFreeFunctions::DoFInfo::dof_access_cell &&
          fe_eval.is_interior_face() == false)
        fe_face_evaluation_process_and_io<VectorizedArrayType::size()>(
          p, n_components, integration_flag, dst_ptr, sm_ptr, fe_eval, temp);
      else
        fe_face_evaluation_process_and_io<1>(
          p, n_components, integration_flag, dst_ptr, sm_ptr, fe_eval, temp);

      return false;
    }

  private:
    template <int fe_degree>
    struct Processor
    {
      static const bool do_integrate = true;
      static const int  dim_         = dim;
      static const int  fe_degree_   = fe_degree;
      using VectorizedArrayType_     = VectorizedArrayType;
      using Number_                  = Number;
      using Number2_                 = Number2;

      template <typename T0, typename T1, typename T2, typename T3, typename T4>
      void
      hermite_grad_vectorized(const T0 &temp_1,
                              const T1 &temp_2,
                              T2        dst_ptr_1,
                              T3        dst_ptr_2,
                              const T4 &grad_weight)
      {
        // case 1a)
        const VectorizedArrayType val  = temp_1 - grad_weight * temp_2;
        const VectorizedArrayType grad = grad_weight * temp_2;
        do_vectorized_add(val, dst_ptr_1);
        do_vectorized_add(grad, dst_ptr_2);
      }

      template <typename T0, typename T1>
      void
      value_vectorized(const T0 &temp, T1 dst_ptr)
      {
        // case 1b)
        do_vectorized_add(temp, dst_ptr);
      }

      template <typename T0, typename T1, typename T2, typename T3>
      void
      hermite_grad_vectorized_indexed(const T0 &temp_1,
                                      const T0 &temp_2,
                                      T1        dst_ptr_1,
                                      T1        dst_ptr_2,
                                      const T2 &grad_weight,
                                      const T3 &indices_1,
                                      const T3 &indices_2)
      {
        // case 2a)
        const VectorizedArrayType val  = temp_1 - grad_weight * temp_2;
        const VectorizedArrayType grad = grad_weight * temp_2;
        do_vectorized_scatter_add(val, indices_1, dst_ptr_1);
        do_vectorized_scatter_add(grad, indices_2, dst_ptr_2);
      }

      template <typename T0, typename T1, typename T2>
      void
      value_vectorized_indexed(const T0 &temp, T1 dst_ptr, const T2 &indices)
      {
        // case 2b)
        do_vectorized_scatter_add(temp, indices, dst_ptr);
      }

      template <typename T0, typename T1, typename T2>
      void
      hermite_grad(const T0 &temp_1,
                   const T0 &temp_2,
                   T1 &      dst_ptr_1,
                   T1 &      dst_ptr_2,
                   const T2 &grad_weight)
      {
        // case 3a)
        const Number val  = temp_1 - grad_weight * temp_2;
        const Number grad = grad_weight * temp_2;
        dst_ptr_1 += val;
        dst_ptr_2 += grad;
      }

      template <typename T0, typename T1>
      void
      value(const T0 &temp, T1 &dst_ptr)
      {
        // case 3b)
        dst_ptr += temp;
      }
    };
  };



  /**
   * This struct implements the action of the inverse mass matrix operation
   * using an FEEvaluationData argument.
   */
  template <int dim, typename Number>
  struct CellwiseInverseMassMatrixImplBasic
  {
    template <int fe_degree, int = 0>
    static bool
    run(const unsigned int                          n_components,
        const FEEvaluationData<dim, Number, false> &fe_eval,
        const Number *                              in_array,
        Number *                                    out_array,
        typename std::enable_if<fe_degree != -1>::type * = nullptr)
    {
      constexpr unsigned int dofs_per_component =
        Utilities::pow(fe_degree + 1, dim);

      Assert(dim >= 1 || dim <= 3, ExcNotImplemented());
      Assert(fe_eval.get_shape_info().element_type <=
               MatrixFreeFunctions::tensor_symmetric,
             ExcNotImplemented());

      EvaluatorTensorProduct<evaluate_evenodd,
                             dim,
                             fe_degree + 1,
                             fe_degree + 1,
                             Number>
        evaluator(
          AlignedVector<Number>(),
          AlignedVector<Number>(),
          fe_eval.get_shape_info().data.front().inverse_shape_values_eo);

      for (unsigned int d = 0; d < n_components; ++d)
        {
          const Number *in  = in_array + d * dofs_per_component;
          Number *      out = out_array + d * dofs_per_component;
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

    template <int fe_degree, int = 0>
    static bool
    run(const unsigned int                          n_components,
        const FEEvaluationData<dim, Number, false> &fe_eval,
        const Number *                              in_array,
        Number *                                    out_array,
        typename std::enable_if<fe_degree == -1>::type * = nullptr)
    {
      static_assert(fe_degree == -1, "Only usable for degree -1");
      const unsigned int dofs_per_component =
        fe_eval.get_shape_info().dofs_per_component_on_cell;

      Assert(dim >= 1 || dim <= 3, ExcNotImplemented());

      EvaluatorTensorProduct<evaluate_general, dim, 0, 0, Number> evaluator(
        fe_eval.get_shape_info().data.front().inverse_shape_values,
        AlignedVector<Number>(),
        AlignedVector<Number>(),
        fe_eval.get_shape_info().data.front().fe_degree + 1,
        fe_eval.get_shape_info().data.front().fe_degree + 1);

      for (unsigned int d = 0; d < n_components; ++d)
        {
          const Number *in  = in_array + d * dofs_per_component;
          Number *      out = out_array + d * dofs_per_component;
          // Need to select 'apply' method with hessian slot because values
          // assume symmetries that do not exist in the inverse shapes
          evaluator.template values<0, true, false>(in, out);
          if (dim > 1)
            evaluator.template values<1, true, false>(out, out);
          if (dim > 2)
            evaluator.template values<2, true, false>(out, out);
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
            evaluator.template values<2, false, false>(out, out);
          if (dim > 1)
            evaluator.template values<1, false, false>(out, out);
          evaluator.template values<0, false, false>(out, out);
        }
      return false;
    }
  };



  /**
   * This struct implements the action of the inverse mass matrix operation
   * using an FEEvaluationData argument.
   */
  template <int dim, typename Number>
  struct CellwiseInverseMassMatrixImplFlexible
  {
    template <int fe_degree, int = 0>
    static bool
    run(const unsigned int           n_desired_components,
        const AlignedVector<Number> &inverse_shape,
        const AlignedVector<Number> &inverse_coefficients,
        const Number *               in_array,
        Number *                     out_array,
        typename std::enable_if<fe_degree != -1>::type * = nullptr)
    {
      constexpr unsigned int dofs_per_component =
        Utilities::pow(fe_degree + 1, dim);
      Assert(inverse_coefficients.size() > 0 &&
               inverse_coefficients.size() % dofs_per_component == 0,
             ExcMessage(
               "Expected diagonal to be a multiple of scalar dof per cells"));
      if (inverse_coefficients.size() != dofs_per_component)
        AssertDimension(n_desired_components * dofs_per_component,
                        inverse_coefficients.size());

      Assert(dim >= 1 || dim <= 3, ExcNotImplemented());

      EvaluatorTensorProduct<evaluate_evenodd,
                             dim,
                             fe_degree + 1,
                             fe_degree + 1,
                             Number>
        evaluator(AlignedVector<Number>(),
                  AlignedVector<Number>(),
                  inverse_shape);

      const unsigned int shift_coefficient =
        inverse_coefficients.size() > dofs_per_component ? dofs_per_component :
                                                           0;
      const Number *inv_coefficient = inverse_coefficients.data();
      for (unsigned int d = 0; d < n_desired_components; ++d)
        {
          const Number *in  = in_array + d * dofs_per_component;
          Number *      out = out_array + d * dofs_per_component;
          // Need to select 'apply' method with hessian slot because values
          // assume symmetries that do not exist in the inverse shapes
          evaluator.template hessians<0, true, false>(in, out);
          if (dim > 1)
            evaluator.template hessians<1, true, false>(out, out);
          if (dim > 2)
            evaluator.template hessians<2, true, false>(out, out);

          for (unsigned int q = 0; q < dofs_per_component; ++q)
            out[q] *= inv_coefficient[q];

          if (dim > 2)
            evaluator.template hessians<2, false, false>(out, out);
          if (dim > 1)
            evaluator.template hessians<1, false, false>(out, out);
          evaluator.template hessians<0, false, false>(out, out);

          inv_coefficient += shift_coefficient;
        }
      return false;
    }

    /**
     * Version for degree = -1
     */
    template <int fe_degree, int = 0>
    static bool
    run(const unsigned int,
        const AlignedVector<Number> &,
        const AlignedVector<Number> &,
        const Number *,
        Number *,
        typename std::enable_if<fe_degree == -1>::type * = nullptr)
    {
      static_assert(fe_degree == -1, "Only usable for degree -1");
      Assert(false, ExcNotImplemented());
      return false;
    }
  };



  /**
   * This struct implements the action of the inverse mass matrix operation
   * using an FEEvaluationData argument.
   */
  template <int dim, typename Number>
  struct CellwiseInverseMassMatrixImplTransformFromQPoints
  {
    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                          n_desired_components,
        const FEEvaluationData<dim, Number, false> &fe_eval,
        const Number *                              in_array,
        Number *                                    out_array)
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

      EvaluatorTensorProduct<do_inplace ? evaluate_evenodd : evaluate_general,
                             dim,
                             fe_degree + 1,
                             n_q_points_1d,
                             Number>
        evaluator(AlignedVector<Number>(),
                  AlignedVector<Number>(),
                  inverse_shape,
                  fe_eval.get_shape_info().data.front().fe_degree + 1,
                  fe_eval.get_shape_info().data.front().n_q_points_1d);

      for (unsigned int d = 0; d < n_desired_components; ++d)
        {
          const Number *in  = in_array + d * n_q_points;
          Number *      out = out_array + d * dofs_per_component;

          auto temp_1 = do_inplace ? out : fe_eval.get_scratch_data().begin();
          auto temp_2 = do_inplace ?
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
