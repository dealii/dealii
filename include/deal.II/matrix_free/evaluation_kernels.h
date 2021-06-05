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


#ifndef dealii_matrix_free_evaluation_kernels_h
#define dealii_matrix_free_evaluation_kernels_h

#include <deal.II/base/config.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/matrix_free/dof_info.h>
#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>
#include <deal.II/matrix_free/type_traits.h>


DEAL_II_NAMESPACE_OPEN


// forward declaration
template <int, typename, bool, typename>
class FEEvaluationBaseData;



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
    static void
    evaluate(const unsigned int                            n_components,
             const EvaluationFlags::EvaluationFlags        evaluation_flag,
             const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
             const Number *                                values_dofs_actual,
             Number *                                      values_quad,
             Number *                                      gradients_quad,
             Number *                                      hessians_quad,
             Number *                                      scratch_data);

    static void
    integrate(const unsigned int                            n_components,
              const EvaluationFlags::EvaluationFlags        integration_flag,
              const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
              Number *                                      values_dofs_actual,
              Number *                                      values_quad,
              Number *                                      gradients_quad,
              Number *                                      scratch_data,
              const bool add_into_values_array);
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
    evaluate(const unsigned int                            n_components,
             const EvaluationFlags::EvaluationFlags        evaluation_flag,
             const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
             const Number *                                values_dofs_actual,
             Number *                                      values_quad,
             Number *                                      gradients_quad,
             Number *                                      hessians_quad,
             Number *                                      scratch_data);

    static void
    integrate(const unsigned int                            n_components,
              const EvaluationFlags::EvaluationFlags        integration_flag,
              const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
              Number *                                      values_dofs_actual,
              Number *                                      values_quad,
              Number *                                      gradients_quad,
              Number *                                      scratch_data,
              const bool add_into_values_array);
  };



  template <MatrixFreeFunctions::ElementType type,
            int                              dim,
            int                              fe_degree,
            int                              n_q_points_1d,
            typename Number>
  inline void
  FEEvaluationImpl<type, dim, fe_degree, n_q_points_1d, Number>::evaluate(
    const unsigned int                            n_components,
    const EvaluationFlags::EvaluationFlags        evaluation_flag,
    const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
    const Number *                                values_dofs_actual,
    Number *                                      values_quad,
    Number *                                      gradients_quad,
    Number *                                      hessians_quad,
    Number *                                      scratch_data)
  {
    if (evaluation_flag == EvaluationFlags::nothing)
      return;

    const EvaluatorVariant variant =
      EvaluatorSelector<type, (fe_degree + n_q_points_1d > 4)>::variant;
    using Eval = EvaluatorTensorProduct<variant,
                                        dim,
                                        fe_degree + 1,
                                        n_q_points_1d,
                                        Number>;
    Eval eval(variant == evaluate_evenodd ?
                shape_info.data.front().shape_values_eo :
                shape_info.data.front().shape_values,
              variant == evaluate_evenodd ?
                shape_info.data.front().shape_gradients_eo :
                shape_info.data.front().shape_gradients,
              variant == evaluate_evenodd ?
                shape_info.data.front().shape_hessians_eo :
                shape_info.data.front().shape_hessians,
              shape_info.data.front().fe_degree + 1,
              shape_info.data.front().n_q_points_1d);

    const unsigned int temp_size =
      Eval::n_rows_of_product == numbers::invalid_unsigned_int ?
        0 :
        (Eval::n_rows_of_product > Eval::n_columns_of_product ?
           Eval::n_rows_of_product :
           Eval::n_columns_of_product);
    Number *temp1 = scratch_data;
    Number *temp2;
    if (temp_size == 0)
      {
        temp2 = temp1 + std::max(Utilities::fixed_power<dim>(
                                   shape_info.data.front().fe_degree + 1),
                                 Utilities::fixed_power<dim>(
                                   shape_info.data.front().n_q_points_1d));
      }
    else
      {
        temp2 = temp1 + temp_size;
      }

    const unsigned int n_q_points =
      temp_size == 0 ? shape_info.n_q_points : Eval::n_columns_of_product;
    const unsigned int dofs_per_comp =
      (type == MatrixFreeFunctions::truncated_tensor) ?
        Utilities::fixed_power<dim>(shape_info.data.front().fe_degree + 1) :
        shape_info.dofs_per_component_on_cell;
    const Number *values_dofs = values_dofs_actual;
    if (type == MatrixFreeFunctions::truncated_tensor)
      {
        Number *values_dofs_tmp =
          scratch_data + 2 * (std::max(shape_info.dofs_per_component_on_cell,
                                       shape_info.n_q_points));
        const int degree =
          fe_degree != -1 ? fe_degree : shape_info.data.front().fe_degree;
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
                      values_dofs_actual
                        [c * shape_info.dofs_per_component_on_cell + count_p];
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

    switch (dim)
      {
        case 1:
          for (unsigned int c = 0; c < n_components; c++)
            {
              if (evaluation_flag & EvaluationFlags::values)
                eval.template values<0, true, false>(values_dofs, values_quad);
              if (evaluation_flag & EvaluationFlags::gradients)
                eval.template gradients<0, true, false>(values_dofs,
                                                        gradients_quad);
              if (evaluation_flag & EvaluationFlags::hessians)
                eval.template hessians<0, true, false>(values_dofs,
                                                       hessians_quad);

              // advance the next component in 1D array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += n_q_points;
              hessians_quad += n_q_points;
            }
          break;

        case 2:
          for (unsigned int c = 0; c < n_components; c++)
            {
              // grad x
              if (evaluation_flag & EvaluationFlags::gradients)
                {
                  eval.template gradients<0, true, false>(values_dofs, temp1);
                  eval.template values<1, true, false>(temp1, gradients_quad);
                }
              if (evaluation_flag & EvaluationFlags::hessians)
                {
                  // grad xy
                  if (!(evaluation_flag & EvaluationFlags::gradients))
                    eval.template gradients<0, true, false>(values_dofs, temp1);
                  eval.template gradients<1, true, false>(temp1,
                                                          hessians_quad +
                                                            2 * n_q_points);

                  // grad xx
                  eval.template hessians<0, true, false>(values_dofs, temp1);
                  eval.template values<1, true, false>(temp1, hessians_quad);
                }

              // grad y
              eval.template values<0, true, false>(values_dofs, temp1);
              if (evaluation_flag & EvaluationFlags::gradients)
                eval.template gradients<1, true, false>(temp1,
                                                        gradients_quad +
                                                          n_q_points);

              // grad yy
              if (evaluation_flag & EvaluationFlags::hessians)
                eval.template hessians<1, true, false>(temp1,
                                                       hessians_quad +
                                                         n_q_points);

              // val: can use values applied in x
              if (evaluation_flag & EvaluationFlags::values)
                eval.template values<1, true, false>(temp1, values_quad);

              // advance to the next component in 1D array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += 2 * n_q_points;
              hessians_quad += 3 * n_q_points;
            }
          break;

        case 3:
          for (unsigned int c = 0; c < n_components; c++)
            {
              if (evaluation_flag & EvaluationFlags::gradients)
                {
                  // grad x
                  eval.template gradients<0, true, false>(values_dofs, temp1);
                  eval.template values<1, true, false>(temp1, temp2);
                  eval.template values<2, true, false>(temp2, gradients_quad);
                }

              if (evaluation_flag & EvaluationFlags::hessians)
                {
                  // grad xz
                  if (!(evaluation_flag & EvaluationFlags::gradients))
                    {
                      eval.template gradients<0, true, false>(values_dofs,
                                                              temp1);
                      eval.template values<1, true, false>(temp1, temp2);
                    }
                  eval.template gradients<2, true, false>(temp2,
                                                          hessians_quad +
                                                            4 * n_q_points);

                  // grad xy
                  eval.template gradients<1, true, false>(temp1, temp2);
                  eval.template values<2, true, false>(temp2,
                                                       hessians_quad +
                                                         3 * n_q_points);

                  // grad xx
                  eval.template hessians<0, true, false>(values_dofs, temp1);
                  eval.template values<1, true, false>(temp1, temp2);
                  eval.template values<2, true, false>(temp2, hessians_quad);
                }

              // grad y
              eval.template values<0, true, false>(values_dofs, temp1);
              if (evaluation_flag & EvaluationFlags::gradients)
                {
                  eval.template gradients<1, true, false>(temp1, temp2);
                  eval.template values<2, true, false>(temp2,
                                                       gradients_quad +
                                                         n_q_points);
                }

              if (evaluation_flag & EvaluationFlags::hessians)
                {
                  // grad yz
                  if (!(evaluation_flag & EvaluationFlags::gradients))
                    eval.template gradients<1, true, false>(temp1, temp2);
                  eval.template gradients<2, true, false>(temp2,
                                                          hessians_quad +
                                                            5 * n_q_points);

                  // grad yy
                  eval.template hessians<1, true, false>(temp1, temp2);
                  eval.template values<2, true, false>(temp2,
                                                       hessians_quad +
                                                         n_q_points);
                }

              // grad z: can use the values applied in x direction stored in
              // temp1
              eval.template values<1, true, false>(temp1, temp2);
              if (evaluation_flag & EvaluationFlags::gradients)
                eval.template gradients<2, true, false>(temp2,
                                                        gradients_quad +
                                                          2 * n_q_points);

              // grad zz: can use the values applied in x and y direction stored
              // in temp2
              if (evaluation_flag & EvaluationFlags::hessians)
                eval.template hessians<2, true, false>(temp2,
                                                       hessians_quad +
                                                         2 * n_q_points);

              // val: can use the values applied in x & y direction stored in
              // temp2
              if (evaluation_flag & EvaluationFlags::values)
                eval.template values<2, true, false>(temp2, values_quad);

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
        for (unsigned int c = 0; c < n_components; ++c)
          for (unsigned int q = 0; q < shape_info.n_q_points; ++q)
            values_quad[c * shape_info.n_q_points + q] +=
              values_dofs[(c + 1) * shape_info.dofs_per_component_on_cell - 1];
      }
  }



  template <MatrixFreeFunctions::ElementType type,
            int                              dim,
            int                              fe_degree,
            int                              n_q_points_1d,
            typename Number>
  inline void
  FEEvaluationImpl<type, dim, fe_degree, n_q_points_1d, Number>::integrate(
    const unsigned int                            n_components,
    const EvaluationFlags::EvaluationFlags        integration_flag,
    const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
    Number *                                      values_dofs_actual,
    Number *                                      values_quad,
    Number *                                      gradients_quad,
    Number *                                      scratch_data,
    const bool                                    add_into_values_array)
  {
    const EvaluatorVariant variant =
      EvaluatorSelector<type, (fe_degree + n_q_points_1d > 4)>::variant;
    using Eval = EvaluatorTensorProduct<variant,
                                        dim,
                                        fe_degree + 1,
                                        n_q_points_1d,
                                        Number>;
    Eval eval(variant == evaluate_evenodd ?
                shape_info.data.front().shape_values_eo :
                shape_info.data.front().shape_values,
              variant == evaluate_evenodd ?
                shape_info.data.front().shape_gradients_eo :
                shape_info.data.front().shape_gradients,
              variant == evaluate_evenodd ?
                shape_info.data.front().shape_hessians_eo :
                shape_info.data.front().shape_hessians,
              shape_info.data.front().fe_degree + 1,
              shape_info.data.front().n_q_points_1d);

    const unsigned int temp_size =
      Eval::n_rows_of_product == numbers::invalid_unsigned_int ?
        0 :
        (Eval::n_rows_of_product > Eval::n_columns_of_product ?
           Eval::n_rows_of_product :
           Eval::n_columns_of_product);
    Number *temp1 = scratch_data;
    Number *temp2;
    if (temp_size == 0)
      {
        temp2 = temp1 + std::max(Utilities::fixed_power<dim>(
                                   shape_info.data.front().fe_degree + 1),
                                 Utilities::fixed_power<dim>(
                                   shape_info.data.front().n_q_points_1d));
      }
    else
      {
        temp2 = temp1 + temp_size;
      }

    const unsigned int n_q_points =
      temp_size == 0 ? shape_info.n_q_points : Eval::n_columns_of_product;
    const unsigned int dofs_per_comp =
      (type == MatrixFreeFunctions::truncated_tensor) ?
        Utilities::fixed_power<dim>(shape_info.data.front().fe_degree + 1) :
        shape_info.dofs_per_component_on_cell;
    // expand dof_values to tensor product for truncated tensor products
    Number *values_dofs =
      (type == MatrixFreeFunctions::truncated_tensor) ?
        scratch_data + 2 * (std::max(shape_info.dofs_per_component_on_cell,
                                     shape_info.n_q_points)) :
        values_dofs_actual;

    switch (dim)
      {
        case 1:
          for (unsigned int c = 0; c < n_components; c++)
            {
              if (integration_flag & EvaluationFlags::values)
                {
                  if (add_into_values_array == false)
                    eval.template values<0, false, false>(values_quad,
                                                          values_dofs);
                  else
                    eval.template values<0, false, true>(values_quad,
                                                         values_dofs);
                }
              if (integration_flag & EvaluationFlags::gradients)
                {
                  if (integration_flag & EvaluationFlags::values ||
                      add_into_values_array == true)
                    eval.template gradients<0, false, true>(gradients_quad,
                                                            values_dofs);
                  else
                    eval.template gradients<0, false, false>(gradients_quad,
                                                             values_dofs);
                }

              // advance to the next component in 1D array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += n_q_points;
            }
          break;

        case 2:
          for (unsigned int c = 0; c < n_components; c++)
            {
              if ((integration_flag & EvaluationFlags::values) &&
                  !(integration_flag & EvaluationFlags::gradients))
                {
                  eval.template values<1, false, false>(values_quad, temp1);
                  if (add_into_values_array == false)
                    eval.template values<0, false, false>(temp1, values_dofs);
                  else
                    eval.template values<0, false, true>(temp1, values_dofs);
                }
              if (integration_flag & EvaluationFlags::gradients)
                {
                  eval.template gradients<1, false, false>(gradients_quad +
                                                             n_q_points,
                                                           temp1);
                  if (integration_flag & EvaluationFlags::values)
                    eval.template values<1, false, true>(values_quad, temp1);
                  if (add_into_values_array == false)
                    eval.template values<0, false, false>(temp1, values_dofs);
                  else
                    eval.template values<0, false, true>(temp1, values_dofs);
                  eval.template values<1, false, false>(gradients_quad, temp1);
                  eval.template gradients<0, false, true>(temp1, values_dofs);
                }

              // advance to the next component in 1D array
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
              gradients_quad += 2 * n_q_points;
            }
          break;

        case 3:
          for (unsigned int c = 0; c < n_components; c++)
            {
              if ((integration_flag & EvaluationFlags::values) &&
                  !(integration_flag & EvaluationFlags::gradients))
                {
                  eval.template values<2, false, false>(values_quad, temp1);
                  eval.template values<1, false, false>(temp1, temp2);
                  if (add_into_values_array == false)
                    eval.template values<0, false, false>(temp2, values_dofs);
                  else
                    eval.template values<0, false, true>(temp2, values_dofs);
                }
              if (integration_flag & EvaluationFlags::gradients)
                {
                  eval.template gradients<2, false, false>(gradients_quad +
                                                             2 * n_q_points,
                                                           temp1);
                  if (integration_flag & EvaluationFlags::values)
                    eval.template values<2, false, true>(values_quad, temp1);
                  eval.template values<1, false, false>(temp1, temp2);
                  eval.template values<2, false, false>(gradients_quad +
                                                          n_q_points,
                                                        temp1);
                  eval.template gradients<1, false, true>(temp1, temp2);
                  if (add_into_values_array == false)
                    eval.template values<0, false, false>(temp2, values_dofs);
                  else
                    eval.template values<0, false, true>(temp2, values_dofs);
                  eval.template values<2, false, false>(gradients_quad, temp1);
                  eval.template values<1, false, false>(temp1, temp2);
                  eval.template gradients<0, false, true>(temp2, values_dofs);
                }

              // advance to the next component in 1D array
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
        values_dofs -= n_components * dofs_per_comp -
                       shape_info.dofs_per_component_on_cell + 1;
        values_quad -= n_components * n_q_points;
        if (integration_flag & EvaluationFlags::values)
          for (unsigned int c = 0; c < n_components; ++c)
            {
              values_dofs[0] = values_quad[0];
              for (unsigned int q = 1; q < shape_info.n_q_points; ++q)
                values_dofs[0] += values_quad[q];
              values_dofs += dofs_per_comp;
              values_quad += n_q_points;
            }
        else
          {
            for (unsigned int c = 0; c < n_components; ++c)
              values_dofs[c * shape_info.dofs_per_component_on_cell] = Number();
            values_dofs += n_components * shape_info.dofs_per_component_on_cell;
          }
      }

    if (type == MatrixFreeFunctions::truncated_tensor)
      {
        values_dofs -= dofs_per_comp * n_components;
        const int degree =
          fe_degree != -1 ? fe_degree : shape_info.data.front().fe_degree;
        for (unsigned int c = 0; c < n_components; ++c)
          for (int i = 0, count_p = 0, count_q = 0;
               i < (dim > 2 ? degree + 1 : 1);
               ++i)
            {
              for (int j = 0; j < (dim > 1 ? degree + 1 - i : 1); ++j)
                {
                  for (int k = 0; k < degree + 1 - j - i;
                       ++k, ++count_p, ++count_q)
                    values_dofs_actual[c *
                                         shape_info.dofs_per_component_on_cell +
                                       count_p] =
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
                      const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                      const Number *values_dofs_actual,
                      Number *      values_quad,
                      Number *      gradients_quad,
                      Number *      hessians_quad,
                      Number *      scratch_data)
  {
    (void)scratch_data;

    const unsigned int n_dofs     = shape_info.dofs_per_component_on_cell;
    const unsigned int n_q_points = shape_info.n_q_points;

    using Eval =
      EvaluatorTensorProduct<evaluate_general, 1, 0, 0, Number, Number>;

    if (evaluation_flag & EvaluationFlags::values)
      {
        const auto shape_values = shape_info.data.front().shape_values.data();
        auto       values_quad_ptr        = values_quad;
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
        const auto shape_gradients =
          shape_info.data.front().shape_gradients.data();
        auto gradients_quad_ptr     = gradients_quad;
        auto values_dofs_actual_ptr = values_dofs_actual;

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
      {
        Assert(false, ExcNotImplemented());
        (void)hessians_quad;
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
                       const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                       Number *   values_dofs_actual,
                       Number *   values_quad,
                       Number *   gradients_quad,
                       Number *   scratch_data,
                       const bool add_into_values_array)
  {
    (void)scratch_data;

    const unsigned int n_dofs     = shape_info.dofs_per_component_on_cell;
    const unsigned int n_q_points = shape_info.n_q_points;

    using Eval =
      EvaluatorTensorProduct<evaluate_general, 1, 0, 0, Number, Number>;

    if (integration_flag & EvaluationFlags::values)
      {
        const auto shape_values = shape_info.data.front().shape_values.data();
        auto       values_quad_ptr        = values_quad;
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
        const auto shape_gradients =
          shape_info.data.front().shape_gradients.data();
        auto gradients_quad_ptr     = gradients_quad;
        auto values_dofs_actual_ptr = values_dofs_actual;

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
                     (integration_flag & EvaluationFlags::values) == false) &&
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
    evaluate(const unsigned int                            n_components,
             const EvaluationFlags::EvaluationFlags        evaluation_flag,
             const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
             const Number *                                values_dofs,
             Number *                                      values_quad,
             Number *                                      gradients_quad,
             Number *                                      hessians_quad,
             Number *                                      scratch_data);

    static void
    integrate(const unsigned int                            n_components,
              const EvaluationFlags::EvaluationFlags        integration_flag,
              const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
              Number *                                      values_dofs,
              Number *                                      values_quad,
              Number *                                      gradients_quad,
              Number *                                      scratch_data,
              const bool add_into_values_array);
  };



  template <int dim, int fe_degree, typename Number>
  inline void
  FEEvaluationImplCollocation<dim, fe_degree, Number>::evaluate(
    const unsigned int                            n_components,
    const EvaluationFlags::EvaluationFlags        evaluation_flag,
    const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
    const Number *                                values_dofs,
    Number *                                      values_quad,
    Number *                                      gradients_quad,
    Number *                                      hessians_quad,
    Number *)
  {
    AssertDimension(
      shape_info.data.front().shape_gradients_collocation_eo.size(),
      (fe_degree + 2) / 2 * (fe_degree + 1));

    EvaluatorTensorProduct<evaluate_evenodd,
                           dim,
                           fe_degree + 1,
                           fe_degree + 1,
                           Number>
                           eval(AlignedVector<Number>(),
           shape_info.data.front().shape_gradients_collocation_eo,
           shape_info.data.front().shape_hessians_collocation_eo);
    constexpr unsigned int n_q_points = Utilities::pow(fe_degree + 1, dim);

    for (unsigned int c = 0; c < n_components; c++)
      {
        if (evaluation_flag & EvaluationFlags::values)
          for (unsigned int i = 0; i < n_q_points; ++i)
            values_quad[i] = values_dofs[i];
        if (evaluation_flag &
            (EvaluationFlags::gradients | EvaluationFlags::hessians))
          {
            eval.template gradients<0, true, false>(values_dofs,
                                                    gradients_quad);
            if (dim > 1)
              eval.template gradients<1, true, false>(values_dofs,
                                                      gradients_quad +
                                                        n_q_points);
            if (dim > 2)
              eval.template gradients<2, true, false>(values_dofs,
                                                      gradients_quad +
                                                        2 * n_q_points);
          }
        if (evaluation_flag & EvaluationFlags::hessians)
          {
            eval.template hessians<0, true, false>(values_dofs, hessians_quad);
            if (dim > 1)
              {
                eval.template gradients<1, true, false>(gradients_quad,
                                                        hessians_quad +
                                                          dim * n_q_points);
                eval.template hessians<1, true, false>(values_dofs,
                                                       hessians_quad +
                                                         n_q_points);
              }
            if (dim > 2)
              {
                eval.template gradients<2, true, false>(gradients_quad,
                                                        hessians_quad +
                                                          4 * n_q_points);
                eval.template gradients<2, true, false>(
                  gradients_quad + n_q_points, hessians_quad + 5 * n_q_points);
                eval.template hessians<2, true, false>(values_dofs,
                                                       hessians_quad +
                                                         2 * n_q_points);
              }
            hessians_quad += (dim * (dim + 1)) / 2 * n_q_points;
          }
        gradients_quad += dim * n_q_points;
        values_quad += n_q_points;
        values_dofs += n_q_points;
      }
  }



  template <int dim, int fe_degree, typename Number>
  inline void
  FEEvaluationImplCollocation<dim, fe_degree, Number>::integrate(
    const unsigned int                            n_components,
    const EvaluationFlags::EvaluationFlags        integration_flag,
    const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
    Number *                                      values_dofs,
    Number *                                      values_quad,
    Number *                                      gradients_quad,
    Number *,
    const bool add_into_values_array)
  {
    AssertDimension(
      shape_info.data.front().shape_gradients_collocation_eo.size(),
      (fe_degree + 2) / 2 * (fe_degree + 1));

    EvaluatorTensorProduct<evaluate_evenodd,
                           dim,
                           fe_degree + 1,
                           fe_degree + 1,
                           Number>
                           eval(AlignedVector<Number>(),
           shape_info.data.front().shape_gradients_collocation_eo,
           shape_info.data.front().shape_hessians_collocation_eo);
    constexpr unsigned int n_q_points = Utilities::pow(fe_degree + 1, dim);

    for (unsigned int c = 0; c < n_components; c++)
      {
        if (integration_flag & EvaluationFlags::values)
          {
            if (add_into_values_array == false)
              for (unsigned int i = 0; i < n_q_points; ++i)
                values_dofs[i] = values_quad[i];
            else
              for (unsigned int i = 0; i < n_q_points; ++i)
                values_dofs[i] += values_quad[i];
          }
        if (integration_flag & EvaluationFlags::gradients)
          {
            if (integration_flag & EvaluationFlags::values ||
                add_into_values_array == true)
              eval.template gradients<0, false, true>(gradients_quad,
                                                      values_dofs);
            else
              eval.template gradients<0, false, false>(gradients_quad,
                                                       values_dofs);
            if (dim > 1)
              eval.template gradients<1, false, true>(gradients_quad +
                                                        n_q_points,
                                                      values_dofs);
            if (dim > 2)
              eval.template gradients<2, false, true>(gradients_quad +
                                                        2 * n_q_points,
                                                      values_dofs);
          }
        gradients_quad += dim * n_q_points;
        values_quad += n_q_points;
        values_dofs += n_q_points;
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
    evaluate(const unsigned int                            n_components,
             const EvaluationFlags::EvaluationFlags        evaluation_flag,
             const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
             const Number *                                values_dofs,
             Number *                                      values_quad,
             Number *                                      gradients_quad,
             Number *                                      hessians_quad,
             Number *                                      scratch_data);

    static void
    integrate(const unsigned int                            n_components,
              const EvaluationFlags::EvaluationFlags        evaluation_flag,
              const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
              Number *                                      values_dofs,
              Number *                                      values_quad,
              Number *                                      gradients_quad,
              Number *                                      scratch_data,
              const bool add_into_values_array);
  };



  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  inline void
  FEEvaluationImplTransformToCollocation<
    dim,
    fe_degree,
    n_q_points_1d,
    Number>::evaluate(const unsigned int                     n_components,
                      const EvaluationFlags::EvaluationFlags evaluation_flag,
                      const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                      const Number *                                values_dofs,
                      Number *                                      values_quad,
                      Number *gradients_quad,
                      Number *hessians_quad,
                      Number *)
  {
    Assert(n_q_points_1d > fe_degree,
           ExcMessage("You lose information when going to a collocation space "
                      "of lower degree, so the evaluation results would be "
                      "wrong. Thus, this class does not permit the desired "
                      "operation."));
    constexpr unsigned int n_q_points = Utilities::pow(n_q_points_1d, dim);

    for (unsigned int c = 0; c < n_components; c++)
      {
        FEEvaluationImplBasisChange<
          evaluate_evenodd,
          EvaluatorQuantity::value,
          dim,
          (fe_degree >= n_q_points_1d ? n_q_points_1d : fe_degree + 1),
          n_q_points_1d,
          Number,
          Number>::do_forward(1,
                              shape_info.data.front().shape_values_eo,
                              values_dofs,
                              values_quad);

        // apply derivatives in the collocation space
        if (evaluation_flag &
            (EvaluationFlags::gradients | EvaluationFlags::hessians))
          FEEvaluationImplCollocation<dim, n_q_points_1d - 1, Number>::evaluate(
            1,
            evaluation_flag &
              (EvaluationFlags::gradients | EvaluationFlags::hessians),
            shape_info,
            values_quad,
            nullptr,
            gradients_quad,
            hessians_quad,
            nullptr);

        values_dofs += shape_info.dofs_per_component_on_cell;
        values_quad += n_q_points;
        gradients_quad += dim * n_q_points;
        hessians_quad += (dim * (dim + 1)) / 2 * n_q_points;
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
                       const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                       Number *values_dofs,
                       Number *values_quad,
                       Number *gradients_quad,
                       Number *,
                       const bool add_into_values_array)
  {
    Assert(n_q_points_1d > fe_degree,
           ExcMessage("You lose information when going to a collocation space "
                      "of lower degree, so the evaluation results would be "
                      "wrong. Thus, this class does not permit the desired "
                      "operation."));
    AssertDimension(
      shape_info.data.front().shape_gradients_collocation_eo.size(),
      (n_q_points_1d + 1) / 2 * n_q_points_1d);
    constexpr unsigned int n_q_points = Utilities::pow(n_q_points_1d, dim);

    for (unsigned int c = 0; c < n_components; c++)
      {
        // apply derivatives in collocation space
        if (integration_flag & EvaluationFlags::gradients)
          FEEvaluationImplCollocation<dim, n_q_points_1d - 1, Number>::
            integrate(1,
                      integration_flag & EvaluationFlags::gradients,
                      shape_info,
                      values_quad,
                      nullptr,
                      gradients_quad,
                      nullptr,
                      /*add_into_values_array=*/integration_flag &
                        EvaluationFlags::values);

        // transform back to the original space
        FEEvaluationImplBasisChange<
          evaluate_evenodd,
          EvaluatorQuantity::value,
          dim,
          (fe_degree >= n_q_points_1d ? n_q_points_1d : fe_degree + 1),
          n_q_points_1d,
          Number,
          Number>::do_backward(1,
                               shape_info.data.front().shape_values_eo,
                               add_into_values_array,
                               values_quad,
                               values_dofs);
        gradients_quad += dim * n_q_points;
        values_quad += n_q_points;
        values_dofs += shape_info.dofs_per_component_on_cell;
      }
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
    run(const unsigned int                                      n_components,
        const EvaluationFlags::EvaluationFlags                  evaluation_flag,
        const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
        Number *values_dofs_actual,
        Number *values_quad,
        Number *gradients_quad,
        Number *hessians_quad,
        Number *scratch_data)
    {
      // We enable a transformation to collocation for derivatives if it gives
      // correct results (first condition), if it is the most efficient choice
      // in terms of operation counts (second condition) and if we were able to
      // initialize the fields in shape_info.templates.h from the polynomials
      // (third condition).
      static constexpr bool use_collocation =
        n_q_points_1d > fe_degree && n_q_points_1d <= 3 * fe_degree / 2 + 1 &&
        n_q_points_1d < 200;

      if (fe_degree >= 0 && fe_degree + 1 == n_q_points_1d &&
          shape_info.element_type ==
            internal::MatrixFreeFunctions::tensor_symmetric_collocation)
        {
          internal::FEEvaluationImplCollocation<dim, fe_degree, Number>::
            evaluate(n_components,
                     evaluation_flag,
                     shape_info,
                     values_dofs_actual,
                     values_quad,
                     gradients_quad,
                     hessians_quad,
                     scratch_data);
        }
      // '<=' on type means tensor_symmetric or tensor_symmetric_hermite, see
      // shape_info.h for more details
      else if (fe_degree >= 0 && use_collocation &&
               shape_info.element_type <=
                 internal::MatrixFreeFunctions::tensor_symmetric)
        {
          internal::FEEvaluationImplTransformToCollocation<
            dim,
            fe_degree,
            n_q_points_1d,
            Number>::evaluate(n_components,
                              evaluation_flag,
                              shape_info,
                              values_dofs_actual,
                              values_quad,
                              gradients_quad,
                              hessians_quad,
                              scratch_data);
        }
      else if (fe_degree >= 0 &&
               shape_info.element_type <=
                 internal::MatrixFreeFunctions::tensor_symmetric)
        {
          internal::FEEvaluationImpl<
            internal::MatrixFreeFunctions::tensor_symmetric,
            dim,
            fe_degree,
            n_q_points_1d,
            Number>::evaluate(n_components,
                              evaluation_flag,
                              shape_info,
                              values_dofs_actual,
                              values_quad,
                              gradients_quad,
                              hessians_quad,
                              scratch_data);
        }
      else if (shape_info.element_type ==
               internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0)
        {
          internal::FEEvaluationImpl<
            internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0,
            dim,
            fe_degree,
            n_q_points_1d,
            Number>::evaluate(n_components,
                              evaluation_flag,
                              shape_info,
                              values_dofs_actual,
                              values_quad,
                              gradients_quad,
                              hessians_quad,
                              scratch_data);
        }
      else if (shape_info.element_type ==
               internal::MatrixFreeFunctions::truncated_tensor)
        {
          internal::FEEvaluationImpl<
            internal::MatrixFreeFunctions::truncated_tensor,
            dim,
            fe_degree,
            n_q_points_1d,
            Number>::evaluate(n_components,
                              evaluation_flag,
                              shape_info,
                              values_dofs_actual,
                              values_quad,
                              gradients_quad,
                              hessians_quad,
                              scratch_data);
        }
      else if (shape_info.element_type ==
               internal::MatrixFreeFunctions::tensor_none)
        {
          internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_none,
                                     dim,
                                     fe_degree,
                                     n_q_points_1d,
                                     Number>::evaluate(n_components,
                                                       evaluation_flag,
                                                       shape_info,
                                                       values_dofs_actual,
                                                       values_quad,
                                                       gradients_quad,
                                                       hessians_quad,
                                                       scratch_data);
        }
      else
        {
          internal::FEEvaluationImpl<
            internal::MatrixFreeFunctions::tensor_general,
            dim,
            fe_degree,
            n_q_points_1d,
            Number>::evaluate(n_components,
                              evaluation_flag,
                              shape_info,
                              values_dofs_actual,
                              values_quad,
                              gradients_quad,
                              hessians_quad,
                              scratch_data);
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
        const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
        Number *   values_dofs_actual,
        Number *   values_quad,
        Number *   gradients_quad,
        Number *   scratch_data,
        const bool sum_into_values_array)
    {
      // We enable a transformation to collocation for derivatives if it gives
      // correct results (first condition), if it is the most efficient choice
      // in terms of operation counts (second condition) and if we were able to
      // initialize the fields in shape_info.templates.h from the polynomials
      // (third condition).
      constexpr bool use_collocation = n_q_points_1d > fe_degree &&
                                       n_q_points_1d <= 3 * fe_degree / 2 + 1 &&
                                       n_q_points_1d < 200;

      if (fe_degree >= 0 && fe_degree + 1 == n_q_points_1d &&
          shape_info.element_type ==
            internal::MatrixFreeFunctions::tensor_symmetric_collocation)
        {
          internal::FEEvaluationImplCollocation<dim, fe_degree, Number>::
            integrate(n_components,
                      integration_flag,
                      shape_info,
                      values_dofs_actual,
                      values_quad,
                      gradients_quad,
                      scratch_data,
                      sum_into_values_array);
        }
      // '<=' on type means tensor_symmetric or tensor_symmetric_hermite, see
      // shape_info.h for more details
      else if (fe_degree >= 0 && use_collocation &&
               shape_info.element_type <=
                 internal::MatrixFreeFunctions::tensor_symmetric)
        {
          internal::FEEvaluationImplTransformToCollocation<
            dim,
            fe_degree,
            n_q_points_1d,
            Number>::integrate(n_components,
                               integration_flag,
                               shape_info,
                               values_dofs_actual,
                               values_quad,
                               gradients_quad,
                               scratch_data,
                               sum_into_values_array);
        }
      else if (fe_degree >= 0 &&
               shape_info.element_type <=
                 internal::MatrixFreeFunctions::tensor_symmetric)
        {
          internal::FEEvaluationImpl<
            internal::MatrixFreeFunctions::tensor_symmetric,
            dim,
            fe_degree,
            n_q_points_1d,
            Number>::integrate(n_components,
                               integration_flag,
                               shape_info,
                               values_dofs_actual,
                               values_quad,
                               gradients_quad,
                               scratch_data,
                               sum_into_values_array);
        }
      else if (shape_info.element_type ==
               internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0)
        {
          internal::FEEvaluationImpl<
            internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0,
            dim,
            fe_degree,
            n_q_points_1d,
            Number>::integrate(n_components,
                               integration_flag,
                               shape_info,
                               values_dofs_actual,
                               values_quad,
                               gradients_quad,
                               scratch_data,
                               sum_into_values_array);
        }
      else if (shape_info.element_type ==
               internal::MatrixFreeFunctions::truncated_tensor)
        {
          internal::FEEvaluationImpl<
            internal::MatrixFreeFunctions::truncated_tensor,
            dim,
            fe_degree,
            n_q_points_1d,
            Number>::integrate(n_components,
                               integration_flag,
                               shape_info,
                               values_dofs_actual,
                               values_quad,
                               gradients_quad,
                               scratch_data,
                               sum_into_values_array);
        }
      else if (shape_info.element_type ==
               internal::MatrixFreeFunctions::tensor_none)
        {
          internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_none,
                                     dim,
                                     fe_degree,
                                     n_q_points_1d,
                                     Number>::integrate(n_components,
                                                        integration_flag,
                                                        shape_info,
                                                        values_dofs_actual,
                                                        values_quad,
                                                        gradients_quad,
                                                        scratch_data,
                                                        sum_into_values_array);
        }
      else
        {
          internal::FEEvaluationImpl<
            internal::MatrixFreeFunctions::tensor_general,
            dim,
            fe_degree,
            n_q_points_1d,
            Number>::integrate(n_components,
                               integration_flag,
                               shape_info,
                               values_dofs_actual,
                               values_quad,
                               gradients_quad,
                               scratch_data,
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
    static constexpr bool use_collocation =
      symmetric_evaluate &&
      n_q_points_1d > fe_degree &&n_q_points_1d <= 3 * fe_degree / 2 + 1 &&
      n_q_points_1d < 200;

    static void
    evaluate_in_face(const unsigned int                            n_components,
                     const MatrixFreeFunctions::ShapeInfo<Number> &data,
                     Number *                                      values_dofs,
                     Number *                                      values_quad,
                     Number *           gradients_quad,
                     Number *           scratch_data,
                     const bool         evaluate_val,
                     const bool         evaluate_grad,
                     const unsigned int subface_index)
    {
      const AlignedVector<Number> &val1 =
        symmetric_evaluate ?
          data.data.front().shape_values_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
             data.data.front().shape_values :
             data.data.front().values_within_subface[subface_index % 2]);
      const AlignedVector<Number> &val2 =
        symmetric_evaluate ?
          data.data.front().shape_values_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
             data.data.front().shape_values :
             data.data.front().values_within_subface[subface_index / 2]);

      const AlignedVector<Number> &grad1 =
        symmetric_evaluate ?
          data.data.front().shape_gradients_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
             data.data.front().shape_gradients :
             data.data.front().gradients_within_subface[subface_index % 2]);
      const AlignedVector<Number> &grad2 =
        symmetric_evaluate ?
          data.data.front().shape_gradients_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
             data.data.front().shape_gradients :
             data.data.front().gradients_within_subface[subface_index / 2]);

      using Eval =
        internal::EvaluatorTensorProduct<symmetric_evaluate ?
                                           internal::evaluate_evenodd :
                                           internal::evaluate_general,
                                         dim - 1,
                                         fe_degree + 1,
                                         n_q_points_1d,
                                         Number>;
      Eval eval1(val1,
                 grad1,
                 AlignedVector<Number>(),
                 data.data.front().fe_degree + 1,
                 data.data.front().n_q_points_1d);
      Eval eval2(val2,
                 grad2,
                 AlignedVector<Number>(),
                 data.data.front().fe_degree + 1,
                 data.data.front().n_q_points_1d);

      const unsigned int size_deg =
        fe_degree > -1 ?
          Utilities::pow(fe_degree + 1, dim - 1) :
          (dim > 1 ?
             Utilities::fixed_power<dim - 1>(data.data.front().fe_degree + 1) :
             1);

      const unsigned int n_q_points = fe_degree > -1 ?
                                        Utilities::pow(n_q_points_1d, dim - 1) :
                                        data.n_q_points_face;

      if (evaluate_grad == false)
        for (unsigned int c = 0; c < n_components; ++c)
          {
            switch (dim)
              {
                case 3:
                  eval1.template values<0, true, false>(values_dofs,
                                                        values_quad);
                  eval2.template values<1, true, false>(values_quad,
                                                        values_quad);
                  break;
                case 2:
                  eval1.template values<0, true, false>(values_dofs,
                                                        values_quad);
                  break;
                case 1:
                  values_quad[0] = values_dofs[0];
                  break;
                default:
                  Assert(false, ExcNotImplemented());
              }
            values_dofs += 2 * size_deg;
            values_quad += n_q_points;
          }
      else
        for (unsigned int c = 0; c < n_components; ++c)
          {
            switch (dim)
              {
                case 3:
                  if (use_collocation)
                    {
                      eval1.template values<0, true, false>(values_dofs,
                                                            values_quad);
                      eval1.template values<1, true, false>(values_quad,
                                                            values_quad);
                      internal::EvaluatorTensorProduct<
                        internal::evaluate_evenodd,
                        dim - 1,
                        n_q_points_1d,
                        n_q_points_1d,
                        Number>
                        eval_grad(
                          AlignedVector<Number>(),
                          data.data.front().shape_gradients_collocation_eo,
                          AlignedVector<Number>());
                      eval_grad.template gradients<0, true, false>(
                        values_quad, gradients_quad);
                      eval_grad.template gradients<1, true, false>(
                        values_quad, gradients_quad + n_q_points);
                    }
                  else
                    {
                      eval1.template gradients<0, true, false>(values_dofs,
                                                               scratch_data);
                      eval2.template values<1, true, false>(scratch_data,
                                                            gradients_quad);

                      eval1.template values<0, true, false>(values_dofs,
                                                            scratch_data);
                      eval2.template gradients<1, true, false>(scratch_data,
                                                               gradients_quad +
                                                                 n_q_points);

                      if (evaluate_val == true)
                        eval2.template values<1, true, false>(scratch_data,
                                                              values_quad);
                    }
                  eval1.template values<0, true, false>(values_dofs + size_deg,
                                                        scratch_data);
                  eval2.template values<1, true, false>(
                    scratch_data, gradients_quad + (dim - 1) * n_q_points);

                  break;
                case 2:
                  eval1.template values<0, true, false>(values_dofs + size_deg,
                                                        gradients_quad +
                                                          (dim - 1) *
                                                            n_q_points);
                  eval1.template gradients<0, true, false>(values_dofs,
                                                           gradients_quad);
                  if (evaluate_val == true)
                    eval1.template values<0, true, false>(values_dofs,
                                                          values_quad);
                  break;
                case 1:
                  values_quad[0]    = values_dofs[0];
                  gradients_quad[0] = values_dofs[1];
                  break;
                default:
                  AssertThrow(false, ExcNotImplemented());
              }
            values_dofs += 2 * size_deg;
            values_quad += n_q_points;
            gradients_quad += dim * n_q_points;
          }
    }

    static void
    integrate_in_face(const unsigned int n_components,
                      const MatrixFreeFunctions::ShapeInfo<Number> &data,
                      Number *                                      values_dofs,
                      Number *                                      values_quad,
                      Number *           gradients_quad,
                      Number *           scratch_data,
                      const bool         integrate_val,
                      const bool         integrate_grad,
                      const unsigned int subface_index)
    {
      const AlignedVector<Number> &val1 =
        symmetric_evaluate ?
          data.data.front().shape_values_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
             data.data.front().shape_values :
             data.data.front().values_within_subface[subface_index % 2]);
      const AlignedVector<Number> &val2 =
        symmetric_evaluate ?
          data.data.front().shape_values_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
             data.data.front().shape_values :
             data.data.front().values_within_subface[subface_index / 2]);

      const AlignedVector<Number> &grad1 =
        symmetric_evaluate ?
          data.data.front().shape_gradients_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
             data.data.front().shape_gradients :
             data.data.front().gradients_within_subface[subface_index % 2]);
      const AlignedVector<Number> &grad2 =
        symmetric_evaluate ?
          data.data.front().shape_gradients_eo :
          (subface_index >= GeometryInfo<dim>::max_children_per_cell ?
             data.data.front().shape_gradients :
             data.data.front().gradients_within_subface[subface_index / 2]);

      using Eval =
        internal::EvaluatorTensorProduct<symmetric_evaluate ?
                                           internal::evaluate_evenodd :
                                           internal::evaluate_general,
                                         dim - 1,
                                         fe_degree + 1,
                                         n_q_points_1d,
                                         Number>;
      Eval eval1(val1,
                 grad1,
                 val1,
                 data.data.front().fe_degree + 1,
                 data.data.front().n_q_points_1d);
      Eval eval2(val2,
                 grad2,
                 val1,
                 data.data.front().fe_degree + 1,
                 data.data.front().n_q_points_1d);

      const unsigned int size_deg =
        fe_degree > -1 ?
          Utilities::pow(fe_degree + 1, dim - 1) :
          (dim > 1 ?
             Utilities::fixed_power<dim - 1>(data.data.front().fe_degree + 1) :
             1);

      const unsigned int n_q_points = fe_degree > -1 ?
                                        Utilities::pow(n_q_points_1d, dim - 1) :
                                        data.n_q_points_face;

      if (integrate_grad == false)
        for (unsigned int c = 0; c < n_components; ++c)
          {
            switch (dim)
              {
                case 3:
                  eval2.template values<1, false, false>(values_quad,
                                                         values_quad);
                  eval1.template values<0, false, false>(values_quad,
                                                         values_dofs);
                  break;
                case 2:
                  eval1.template values<0, false, false>(values_quad,
                                                         values_dofs);
                  break;
                case 1:
                  values_dofs[0] = values_quad[0];
                  break;
                default:
                  Assert(false, ExcNotImplemented());
              }
            values_dofs += 2 * size_deg;
            values_quad += n_q_points;
          }
      else
        for (unsigned int c = 0; c < n_components; ++c)
          {
            switch (dim)
              {
                case 3:
                  eval2.template values<1, false, false>(gradients_quad +
                                                           2 * n_q_points,
                                                         gradients_quad +
                                                           2 * n_q_points);
                  eval1.template values<0, false, false>(
                    gradients_quad + 2 * n_q_points, values_dofs + size_deg);
                  if (use_collocation)
                    {
                      internal::EvaluatorTensorProduct<
                        internal::evaluate_evenodd,
                        dim - 1,
                        n_q_points_1d,
                        n_q_points_1d,
                        Number>
                        eval_grad(
                          AlignedVector<Number>(),
                          data.data.front().shape_gradients_collocation_eo,
                          AlignedVector<Number>());
                      if (integrate_val)
                        eval_grad.template gradients<1, false, true>(
                          gradients_quad + n_q_points, values_quad);
                      else
                        eval_grad.template gradients<1, false, false>(
                          gradients_quad + n_q_points, values_quad);
                      eval_grad.template gradients<0, false, true>(
                        gradients_quad, values_quad);
                      eval1.template values<1, false, false>(values_quad,
                                                             values_quad);
                      eval1.template values<0, false, false>(values_quad,
                                                             values_dofs);
                    }
                  else
                    {
                      if (integrate_val)
                        {
                          eval2.template values<1, false, false>(values_quad,
                                                                 scratch_data);
                          eval2.template gradients<1, false, true>(
                            gradients_quad + n_q_points, scratch_data);
                        }
                      else
                        eval2.template gradients<1, false, false>(
                          gradients_quad + n_q_points, scratch_data);

                      eval1.template values<0, false, false>(scratch_data,
                                                             values_dofs);
                      eval2.template values<1, false, false>(gradients_quad,
                                                             scratch_data);
                      eval1.template gradients<0, false, true>(scratch_data,
                                                               values_dofs);
                    }
                  break;
                case 2:
                  eval1.template values<0, false, false>(
                    gradients_quad + n_q_points, values_dofs + size_deg);
                  eval1.template gradients<0, false, false>(gradients_quad,
                                                            values_dofs);
                  if (integrate_val == true)
                    eval1.template values<0, false, true>(values_quad,
                                                          values_dofs);
                  break;
                case 1:
                  values_dofs[0] = values_quad[0];
                  values_dofs[1] = gradients_quad[0];
                  break;
                default:
                  AssertThrow(false, ExcNotImplemented());
              }
            values_dofs += 2 * size_deg;
            values_quad += n_q_points;
            gradients_quad += dim * n_q_points;
          }
    }
  };



  template <int dim, int fe_degree, typename Number, bool lex_faces = false>
  struct FEFaceNormalEvaluationImpl
  {
    template <bool do_evaluate, bool add_into_output>
    static void
    interpolate(const unsigned int                            n_components,
                const MatrixFreeFunctions::ShapeInfo<Number> &data,
                const Number *                                input,
                Number *                                      output,
                const bool                                    do_gradients,
                const unsigned int                            face_no)
    {
      Assert(static_cast<unsigned int>(fe_degree) ==
                 data.data.front().fe_degree ||
               fe_degree == -1,
             ExcInternalError());

      interpolate_generic<do_evaluate, add_into_output>(
        n_components,
        input,
        output,
        do_gradients,
        face_no,
        data.data.front().fe_degree + 1,
        data.data.front().shape_data_on_face,
        data.dofs_per_component_on_cell,
        2 * data.dofs_per_component_on_face);
    }

    /**
     * Interpolate the values on the cell quadrature points onto a face.
     */
    template <bool do_evaluate, bool add_into_output>
    static void
    interpolate_quadrature(const unsigned int n_components,
                           const MatrixFreeFunctions::ShapeInfo<Number> &data,
                           const Number *                                input,
                           Number *                                      output,
                           const bool         do_gradients,
                           const unsigned int face_no)
    {
      Assert(static_cast<unsigned int>(fe_degree + 1) ==
                 data.data.front().quadrature.size() ||
               fe_degree == -1,
             ExcInternalError());

      interpolate_generic<do_evaluate, add_into_output>(
        n_components,
        input,
        output,
        do_gradients,
        face_no,
        data.data.front().quadrature.size(),
        data.data.front().quadrature_data_on_face,
        data.n_q_points,
        data.n_q_points_face);
    }

  private:
    template <bool do_evaluate, bool add_into_output, int face_direction = 0>
    static void
    interpolate_generic(const unsigned int n_components,
                        const Number *     input,
                        Number *           output,
                        const bool         do_gradients,
                        const unsigned int face_no,
                        const unsigned int n_points_1d,
                        const std::array<AlignedVector<Number>, 2> &shape_data,
                        const unsigned int dofs_per_component_on_cell,
                        const unsigned int dofs_per_component_on_face)
    {
      if (face_direction == face_no / 2)
        {
          internal::EvaluatorTensorProduct<internal::evaluate_general,
                                           dim,
                                           fe_degree + 1,
                                           0,
                                           Number>
            evalf(shape_data[face_no % 2],
                  AlignedVector<Number>(),
                  AlignedVector<Number>(),
                  n_points_1d,
                  0);

          const unsigned int in_stride = do_evaluate ?
                                           dofs_per_component_on_cell :
                                           dofs_per_component_on_face;
          const unsigned int out_stride = do_evaluate ?
                                            dofs_per_component_on_face :
                                            dofs_per_component_on_cell;

          for (unsigned int c = 0; c < n_components; c++)
            {
              if (do_gradients)
                evalf.template apply_face<face_direction,
                                          do_evaluate,
                                          add_into_output,
                                          1,
                                          lex_faces>(input, output);
              else
                evalf.template apply_face<face_direction,
                                          do_evaluate,
                                          add_into_output,
                                          0,
                                          lex_faces>(input, output);
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
            do_gradients,
            face_no,
            n_points_1d,
            shape_data,
            dofs_per_component_on_cell,
            dofs_per_component_on_face);
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
  template <typename Number, unsigned int width>
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
  template <typename Number, unsigned int width>
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
  template <typename Number, unsigned int width>
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
  template <typename Number, unsigned int width>
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
  adjust_for_face_orientation(const unsigned int            dim,
                              const unsigned int            n_components,
                              const unsigned int            face_orientation,
                              const Table<2, unsigned int> &orientation_map,
                              const bool                    integrate,
                              const bool                    values,
                              const bool                    gradients,
                              const unsigned int            n_q_points,
                              Number *                      tmp_values,
                              Number *                      values_quad,
                              Number *                      gradients_quad)
  {
    Assert(face_orientation, ExcInternalError());
    const unsigned int *orientation = &orientation_map[face_orientation][0];
    for (unsigned int c = 0; c < n_components; ++c)
      {
        if (values == true)
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
        if (gradients == true)
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
      }
  }



  template <int dim, typename VectorizedArrayType>
  struct FEFaceEvaluationImplEvaluateSelector
  {
    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                                         n_components,
        const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
        const VectorizedArrayType *                                values_array,
        VectorizedArrayType *                                      values_quad,
        VectorizedArrayType *         gradients_quad,
        VectorizedArrayType *         scratch_data,
        const bool                    evaluate_values,
        const bool                    evaluate_gradients,
        const unsigned int            face_no,
        const unsigned int            subface_index,
        const unsigned int            face_orientation,
        const Table<2, unsigned int> &orientation_map)
    {
      if (data.element_type == MatrixFreeFunctions::tensor_none)
        {
          const unsigned int n_dofs     = data.dofs_per_component_on_cell;
          const unsigned int n_q_points = data.n_q_points_faces[face_no];
          const auto         shape_info = data.data.front();

          using Eval = EvaluatorTensorProduct<evaluate_general,
                                              1,
                                              0,
                                              0,
                                              VectorizedArrayType,
                                              VectorizedArrayType>;

          if (evaluate_values)
            {
              const auto shape_values =
                &shape_info.shape_values_face(face_no, face_orientation, 0);

              auto values_quad_ptr        = values_quad;
              auto values_dofs_actual_ptr = values_array;

              Eval eval(shape_values, nullptr, nullptr, n_dofs, n_q_points);
              for (unsigned int c = 0; c < n_components; ++c)
                {
                  eval.template values<0, true, false>(values_dofs_actual_ptr,
                                                       values_quad_ptr);

                  values_quad_ptr += n_q_points;
                  values_dofs_actual_ptr += n_dofs;
                }
            }

          if (evaluate_gradients)
            {
              auto gradients_quad_ptr     = gradients_quad;
              auto values_dofs_actual_ptr = values_array;

              std::array<const VectorizedArrayType *, dim> shape_gradients;
              for (unsigned int d = 0; d < dim; ++d)
                shape_gradients[d] = &shape_info.shape_gradients_face(
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


          return true;
        }

      constexpr unsigned int static_dofs_per_face =
        fe_degree > -1 ? Utilities::pow(fe_degree + 1, dim - 1) :
                         numbers::invalid_unsigned_int;
      const unsigned int dofs_per_face =
        fe_degree > -1 ?
          static_dofs_per_face :
          Utilities::pow(data.data.front().fe_degree + 1, dim - 1);

      VectorizedArrayType *temp1 = scratch_data;

      FEFaceNormalEvaluationImpl<dim, fe_degree, VectorizedArrayType>::
        template interpolate<true, false>(
          n_components, data, values_array, temp1, evaluate_gradients, face_no);

      const unsigned int n_q_points_1d_actual =
        fe_degree > -1 ? n_q_points_1d : 0;
      if (fe_degree > -1 &&
          subface_index >= GeometryInfo<dim>::max_children_per_cell &&
          data.element_type <= MatrixFreeFunctions::tensor_symmetric)
        FEFaceEvaluationImpl<
          true,
          dim,
          fe_degree,
          n_q_points_1d_actual,
          VectorizedArrayType>::evaluate_in_face(n_components,
                                                 data,
                                                 temp1,
                                                 values_quad,
                                                 gradients_quad,
                                                 scratch_data + 2 *
                                                                  n_components *
                                                                  dofs_per_face,
                                                 evaluate_values,
                                                 evaluate_gradients,
                                                 subface_index);
      else
        FEFaceEvaluationImpl<
          false,
          dim,
          fe_degree,
          n_q_points_1d_actual,
          VectorizedArrayType>::evaluate_in_face(n_components,
                                                 data,
                                                 temp1,
                                                 values_quad,
                                                 gradients_quad,
                                                 scratch_data + 2 *
                                                                  n_components *
                                                                  dofs_per_face,
                                                 evaluate_values,
                                                 evaluate_gradients,
                                                 subface_index);

      if (face_orientation)
        adjust_for_face_orientation(dim,
                                    n_components,
                                    face_orientation,
                                    orientation_map,
                                    false,
                                    evaluate_values,
                                    evaluate_gradients,
                                    data.n_q_points_face,
                                    scratch_data,
                                    values_quad,
                                    gradients_quad);

      return false;
    }
  };



  template <int dim, typename VectorizedArrayType>
  struct FEFaceEvaluationImplIntegrateSelector
  {
    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                                         n_components,
        const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
        VectorizedArrayType *                                      values_array,
        VectorizedArrayType *                                      values_quad,
        VectorizedArrayType *         gradients_quad,
        VectorizedArrayType *         scratch_data,
        const bool                    integrate_values,
        const bool                    integrate_gradients,
        const unsigned int            face_no,
        const unsigned int            subface_index,
        const unsigned int            face_orientation,
        const Table<2, unsigned int> &orientation_map)
    {
      if (data.element_type == MatrixFreeFunctions::tensor_none)
        {
          const unsigned int n_dofs     = data.dofs_per_component_on_cell;
          const unsigned int n_q_points = data.n_q_points_faces[face_no];
          const auto         shape_info = data.data.front();

          using Eval = EvaluatorTensorProduct<evaluate_general,
                                              1,
                                              0,
                                              0,
                                              VectorizedArrayType,
                                              VectorizedArrayType>;

          if (integrate_values)
            {
              const auto shape_values =
                &shape_info.shape_values_face(face_no, face_orientation, 0);

              auto values_quad_ptr        = values_quad;
              auto values_dofs_actual_ptr = values_array;

              Eval eval(shape_values, nullptr, nullptr, n_dofs, n_q_points);
              for (unsigned int c = 0; c < n_components; ++c)
                {
                  eval.template values<0, false, false>(values_quad_ptr,
                                                        values_dofs_actual_ptr);

                  values_quad_ptr += n_q_points;
                  values_dofs_actual_ptr += n_dofs;
                }
            }

          if (integrate_gradients)
            {
              auto gradients_quad_ptr     = gradients_quad;
              auto values_dofs_actual_ptr = values_array;

              std::array<const VectorizedArrayType *, dim> shape_gradients;
              for (unsigned int d = 0; d < dim; ++d)
                shape_gradients[d] = &shape_info.shape_gradients_face(
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

                      if ((integrate_values == false) && d == 0)
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


          return true;
        }

      if (face_orientation)
        adjust_for_face_orientation(dim,
                                    n_components,
                                    face_orientation,
                                    orientation_map,
                                    true,
                                    integrate_values,
                                    integrate_gradients,
                                    data.n_q_points_face,
                                    scratch_data,
                                    values_quad,
                                    gradients_quad);

      constexpr unsigned int static_dofs_per_face =
        fe_degree > -1 ? Utilities::pow(fe_degree + 1, dim - 1) :
                         numbers::invalid_unsigned_int;
      const unsigned int dofs_per_face =
        fe_degree > -1 ?
          static_dofs_per_face :
          Utilities::pow(data.data.front().fe_degree + 1, dim - 1);

      VectorizedArrayType *temp1 = scratch_data;

      const unsigned int n_q_points_1d_actual =
        fe_degree > -1 ? n_q_points_1d : 0;
      if (fe_degree > -1 &&
          subface_index >= GeometryInfo<dim - 1>::max_children_per_cell &&
          data.element_type <= MatrixFreeFunctions::tensor_symmetric)
        FEFaceEvaluationImpl<
          true,
          dim,
          fe_degree,
          n_q_points_1d_actual,
          VectorizedArrayType>::integrate_in_face(n_components,
                                                  data,
                                                  temp1,
                                                  values_quad,
                                                  gradients_quad,
                                                  scratch_data +
                                                    2 * n_components *
                                                      dofs_per_face,
                                                  integrate_values,
                                                  integrate_gradients,
                                                  subface_index);
      else
        FEFaceEvaluationImpl<
          false,
          dim,
          fe_degree,
          n_q_points_1d_actual,
          VectorizedArrayType>::integrate_in_face(n_components,
                                                  data,
                                                  temp1,
                                                  values_quad,
                                                  gradients_quad,
                                                  scratch_data +
                                                    2 * n_components *
                                                      dofs_per_face,
                                                  integrate_values,
                                                  integrate_gradients,
                                                  subface_index);

      FEFaceNormalEvaluationImpl<dim, fe_degree, VectorizedArrayType>::
        template interpolate<false, false>(n_components,
                                           data,
                                           temp1,
                                           values_array,
                                           integrate_gradients,
                                           face_no);
      return false;
    }
  };



  template <int n_face_orientations, typename Processor>
  static bool
  fe_face_evaluation_process_and_io(Processor &proc)
  {
    auto  n_components             = proc.n_components;
    auto  integrate                = proc.integrate;
    auto  global_vector_ptr        = proc.global_vector_ptr;
    auto &sm_ptr                   = proc.sm_ptr;
    auto &data                     = proc.data;
    auto &dof_info                 = proc.dof_info;
    auto  values_quad              = proc.values_quad;
    auto  gradients_quad           = proc.gradients_quad;
    auto  scratch_data             = proc.scratch_data;
    auto  do_values                = proc.do_values;
    auto  do_gradients             = proc.do_gradients;
    auto  active_fe_index          = proc.active_fe_index;
    auto  first_selected_component = proc.first_selected_component;
    auto  cells                    = proc.cells;
    auto  face_nos                 = proc.face_nos;
    auto  subface_index            = proc.subface_index;
    auto  dof_access_index         = proc.dof_access_index;
    auto  face_orientations        = proc.face_orientations;
    auto &orientation_map          = proc.orientation_map;

    static const int dim       = Processor::dim_;
    static const int fe_degree = Processor::fe_degree_;
    using VectorizedArrayType  = typename Processor::VectorizedArrayType_;

    using Number   = typename Processor::Number_;
    using Number2_ = typename Processor::Number2_;

    const unsigned int cell = cells[0];

    // In the case of integration, we do not need to reshuffle the
    // data at the quadrature points to adjust for the face
    // orientation if the shape functions are nodal at the cell
    // boundaries (and we only requested the integration of the
    // values) or Hermite shape functions are used. These cases are
    // handled later when the values are written back into the
    // glrobal vector.
    if (integrate &&
        (face_orientations[0] > 0 &&
         (subface_index < GeometryInfo<dim>::max_children_per_cell ||
          !(((do_gradients == false &&
              data.data.front().nodal_at_cell_boundaries == true &&
              fe_degree > 0) ||
             (data.element_type ==
                MatrixFreeFunctions::tensor_symmetric_hermite &&
              fe_degree > 1)) &&
            (dof_info.index_storage_variants[dof_access_index][cell] ==
               MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
                 interleaved_contiguous ||
             dof_info.index_storage_variants[dof_access_index][cell] ==
               MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
                 interleaved_contiguous_strided ||
             dof_info.index_storage_variants[dof_access_index][cell] ==
               MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
                 interleaved_contiguous_mixed_strides ||
             dof_info.index_storage_variants[dof_access_index][cell] ==
               MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
                 contiguous)))))
      {
        AssertDimension(n_face_orientations, 1);
        adjust_for_face_orientation(dim,
                                    n_components,
                                    face_orientations[0],
                                    orientation_map,
                                    true,
                                    do_values,
                                    do_gradients,
                                    data.n_q_points_face,
                                    scratch_data,
                                    values_quad,
                                    gradients_quad);
      }

    // we know that the gradient weights for the Hermite case on the
    // right (side==1) are the negative from the value at the left
    // (side==0), so we only read out one of them.
    VectorizedArrayType grad_weight =
      (data.data.front().nodal_at_cell_boundaries == true && fe_degree > 1 &&
       data.element_type == MatrixFreeFunctions::tensor_symmetric_hermite) ?
        data.data.front()
          .shape_data_on_face[0][fe_degree + (integrate ?
                                                (2 - (face_nos[0] % 2)) :
                                                (1 + (face_nos[0] % 2)))] :
        VectorizedArrayType(0.0 /*dummy*/);

    constexpr unsigned int static_dofs_per_component =
      fe_degree > -1 ? Utilities::pow(fe_degree + 1, dim) :
                       numbers::invalid_unsigned_int;
    constexpr unsigned int static_dofs_per_face =
      fe_degree > -1 ? Utilities::pow(fe_degree + 1, dim - 1) :
                       numbers::invalid_unsigned_int;
    const unsigned int dofs_per_face =
      fe_degree > -1 ? static_dofs_per_face :
                       Utilities::pow(data.data.front().fe_degree + 1, dim - 1);

    VectorizedArrayType *temp1 = scratch_data;

    const unsigned int dummy = 0;

    // re-orientation
    std::array<const unsigned int *, n_face_orientations> orientation = {};

    if (n_face_orientations == 1)
      orientation[0] = (data.data.front().nodal_at_cell_boundaries == true) ?
                         &data.face_orientations[face_orientations[0]][0] :
                         &dummy;
    else
      {
        for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
          {
            // the loop breaks once an invalid_unsigned_int is hit for
            // all cases except the exterior faces in the ECL loop (where
            // some faces might be at the boundaries but others not)
            if (cells[v] == numbers::invalid_unsigned_int)
              continue;

            orientation[v] =
              (data.data.front().nodal_at_cell_boundaries == true) ?
                &data.face_orientations[face_orientations[v]][0] :
                &dummy;
          }
      }

    // face_to_cell_index_hermite
    std::array<const unsigned int *, n_face_orientations> index_array_hermite =
      {};

    if (n_face_orientations == 1)
      index_array_hermite[0] =
        (data.data.front().nodal_at_cell_boundaries == true && fe_degree > 1 &&
         data.element_type == MatrixFreeFunctions::tensor_symmetric_hermite) ?
          &data.face_to_cell_index_hermite(face_nos[0], 0) :
          &dummy;

    if (n_face_orientations > 1 &&
        data.data.front().nodal_at_cell_boundaries == true && fe_degree > 1 &&
        data.element_type == MatrixFreeFunctions::tensor_symmetric_hermite)
      {
        for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
          {
            if (cells[v] == numbers::invalid_unsigned_int)
              continue;

            grad_weight[v] =
              data.data.front().shape_data_on_face
                [0][fe_degree + (integrate ? (2 - (face_nos[v] % 2)) :
                                             (1 + (face_nos[v] % 2)))][v];

            index_array_hermite[v] =
              &data.face_to_cell_index_hermite(face_nos[v], 0);
          }
      }

    // face_to_cell_index_nodal
    std::array<const unsigned int *, n_face_orientations> index_array_nodal =
      {};

    if (n_face_orientations == 1)
      index_array_nodal[0] =
        (data.data.front().nodal_at_cell_boundaries == true) ?
          &data.face_to_cell_index_nodal(face_nos[0], 0) :
          &dummy;

    if (n_face_orientations > 1 &&
        (data.data.front().nodal_at_cell_boundaries == true))
      {
        for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
          {
            if (cells[v] == numbers::invalid_unsigned_int)
              continue;

            index_array_nodal[v] =
              &data.face_to_cell_index_nodal(face_nos[v], 0);
          }
      }

    const auto reorientate = [&](const unsigned int v, const unsigned int i) {
      return (dim < 3 ||
              face_orientations[n_face_orientations == 1 ? 0 : v] == 0 ||
              subface_index < GeometryInfo<dim>::max_children_per_cell) ?
               i :
               orientation[v][i];
    };

    // this variable keeps track of whether we are able to directly write
    // the results into the result (function returns true) or not, requiring
    // an additional call to another function
    bool accesses_global_vector = true;

    for (unsigned int comp = 0; comp < n_components; ++comp)
      {
        if (integrate)
          proc.in_face_operation(temp1, comp);

        // we can only use the fast functions if we know the polynomial degree
        // as a template parameter (fe_degree != -1), and it only makes sense
        // to use the functions for at least linear functions for values on
        // the faces and quadratic functions for gradients on the faces, so
        // include the switch here
        if ((do_gradients == false &&
             data.data.front().nodal_at_cell_boundaries == true &&
             fe_degree > 0) ||
            (data.element_type ==
               MatrixFreeFunctions::tensor_symmetric_hermite &&
             fe_degree > 1))
          {
            // case 1: contiguous and interleaved indices
            if (n_face_orientations == 1 &&
                dof_info.index_storage_variants[dof_access_index][cell] ==
                  MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
                    interleaved_contiguous)
              {
                AssertDimension(n_face_orientations, 1);

                AssertDimension(
                  dof_info.n_vectorization_lanes_filled[dof_access_index][cell],
                  VectorizedArrayType::size());
                Number2_ *vector_ptr =
                  global_vector_ptr +
                  dof_info.dof_indices_contiguous[dof_access_index]
                                                 [cell *
                                                  VectorizedArrayType::size()] +
                  (dof_info
                     .component_dof_indices_offset[active_fe_index]
                                                  [first_selected_component] +
                   comp * static_dofs_per_component) *
                    VectorizedArrayType::size();

                if (fe_degree > 1 && do_gradients == true)
                  {
                    for (unsigned int i = 0; i < dofs_per_face; ++i)
                      {
                        if (n_face_orientations == 1)
                          {
                            const unsigned int ind1 =
                              index_array_hermite[0][2 * i];
                            const unsigned int ind2 =
                              index_array_hermite[0][2 * i + 1];
                            AssertIndexRange(ind1,
                                             data.dofs_per_component_on_cell);
                            AssertIndexRange(ind2,
                                             data.dofs_per_component_on_cell);
                            const unsigned int i_ = reorientate(0, i);
                            proc.hermite_grad_vectorized(
                              temp1[i_],
                              temp1[i_ + dofs_per_face],
                              vector_ptr + ind1 * VectorizedArrayType::size(),
                              vector_ptr + ind2 * VectorizedArrayType::size(),
                              grad_weight);
                          }
                        else
                          {
                            Assert(false, ExcNotImplemented());
                          }
                      }
                  }
                else
                  {
                    for (unsigned int i = 0; i < dofs_per_face; ++i)
                      {
                        if (n_face_orientations == 1)
                          {
                            const unsigned int i_  = reorientate(0, i);
                            const unsigned int ind = index_array_nodal[0][i];
                            proc.value_vectorized(
                              temp1[i_],
                              vector_ptr + ind * VectorizedArrayType::size());
                          }
                        else
                          {
                            Assert(false, ExcNotImplemented());
                          }
                      }
                  }
              }

            // case 2: contiguous and interleaved indices with fixed stride
            else if (n_face_orientations == 1 &&
                     dof_info.index_storage_variants[dof_access_index][cell] ==
                       MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
                         interleaved_contiguous_strided)
              {
                AssertDimension(n_face_orientations, 1);

                AssertDimension(
                  dof_info.n_vectorization_lanes_filled[dof_access_index][cell],
                  VectorizedArrayType::size());
                const unsigned int *indices =
                  &dof_info.dof_indices_contiguous[dof_access_index]
                                                  [cell *
                                                   VectorizedArrayType::size()];
                Number2_ *vector_ptr =
                  global_vector_ptr +
                  (comp * static_dofs_per_component +
                   dof_info
                     .component_dof_indices_offset[active_fe_index]
                                                  [first_selected_component]) *
                    VectorizedArrayType::size();
                if (fe_degree > 1 && do_gradients == true)
                  {
                    for (unsigned int i = 0; i < dofs_per_face; ++i)
                      {
                        if (n_face_orientations == 1)
                          {
                            const unsigned int i_ = reorientate(0, i);
                            const unsigned int ind1 =
                              index_array_hermite[0][2 * i] *
                              VectorizedArrayType::size();
                            const unsigned int ind2 =
                              index_array_hermite[0][2 * i + 1] *
                              VectorizedArrayType::size();
                            proc.hermite_grad_vectorized_indexed(
                              temp1[i_],
                              temp1[i_ + dofs_per_face],
                              vector_ptr + ind1,
                              vector_ptr + ind2,
                              grad_weight,
                              indices,
                              indices);
                          }
                        else
                          {
                            Assert(false, ExcNotImplemented());
                          }
                      }
                  }
                else
                  {
                    for (unsigned int i = 0; i < dofs_per_face; ++i)
                      {
                        if (n_face_orientations == 1)
                          {
                            const unsigned int i_ = reorientate(0, i);
                            const unsigned int ind =
                              index_array_nodal[0][i] *
                              VectorizedArrayType::size();
                            proc.value_vectorized_indexed(temp1[i_],
                                                          vector_ptr + ind,
                                                          indices);
                          }
                        else
                          {
                            Assert(false, ExcNotImplemented());
                          }
                      }
                  }
              }

            // case 3: contiguous and interleaved indices with mixed stride
            else if (n_face_orientations == 1 &&
                     dof_info.index_storage_variants[dof_access_index][cell] ==
                       MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
                         interleaved_contiguous_mixed_strides)
              {
                AssertDimension(n_face_orientations, 1);

                const unsigned int *strides =
                  &dof_info.dof_indices_interleave_strides
                     [dof_access_index][cell * VectorizedArrayType::size()];
                unsigned int indices[VectorizedArrayType::size()];
                for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
                  indices[v] =
                    dof_info.dof_indices_contiguous
                      [dof_access_index]
                      [cell * VectorizedArrayType::size() + v] +
                    (dof_info
                       .component_dof_indices_offset[active_fe_index]
                                                    [first_selected_component] +
                     comp * static_dofs_per_component) *
                      strides[v];
                const unsigned int n_filled_lanes =
                  dof_info.n_vectorization_lanes_filled[dof_access_index][cell];

                if (fe_degree > 1 && do_gradients == true)
                  {
                    if (n_filled_lanes == VectorizedArrayType::size())
                      for (unsigned int i = 0; i < dofs_per_face; ++i)
                        {
                          if (n_face_orientations == 1)
                            {
                              const unsigned int i_ = reorientate(0, i);
                              unsigned int ind1[VectorizedArrayType::size()];
                              DEAL_II_OPENMP_SIMD_PRAGMA
                              for (unsigned int v = 0;
                                   v < VectorizedArrayType::size();
                                   ++v)
                                ind1[v] =
                                  indices[v] +
                                  index_array_hermite[0 /*TODO*/][2 * i] *
                                    strides[v];
                              unsigned int ind2[VectorizedArrayType::size()];
                              DEAL_II_OPENMP_SIMD_PRAGMA
                              for (unsigned int v = 0;
                                   v < VectorizedArrayType::size();
                                   ++v)
                                ind2[v] =
                                  indices[v] +
                                  index_array_hermite[0 /*TODO*/][2 * i + 1] *
                                    strides[v];
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
                              Assert(false, ExcNotImplemented());
                            }
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
                                reorientate(n_face_orientations == 1 ? 0 : v,
                                            i);
                              proc.hermite_grad(
                                temp1[i_][v],
                                temp1[i_ + dofs_per_face][v],
                                global_vector_ptr
                                  [indices[v] +
                                   index_array_hermite
                                       [n_face_orientations == 1 ? 0 : v]
                                       [2 * i] *
                                     strides[v]],
                                global_vector_ptr
                                  [indices[v] +
                                   index_array_hermite
                                       [n_face_orientations == 1 ? 0 : v]
                                       [2 * i + 1] *
                                     strides[v]],
                                grad_weight[n_face_orientations == 1 ? 0 : v]);
                            }
                      }
                  }
                else
                  {
                    if (n_filled_lanes == VectorizedArrayType::size())
                      for (unsigned int i = 0; i < dofs_per_face; ++i)
                        {
                          if (n_face_orientations == 1)
                            {
                              unsigned int ind[VectorizedArrayType::size()];
                              DEAL_II_OPENMP_SIMD_PRAGMA
                              for (unsigned int v = 0;
                                   v < VectorizedArrayType::size();
                                   ++v)
                                ind[v] = indices[v] +
                                         index_array_nodal[0][i] * strides[v];
                              const unsigned int i_ = reorientate(0, i);
                              proc.value_vectorized_indexed(temp1[i_],
                                                            global_vector_ptr,
                                                            ind);
                            }
                          else
                            {
                              Assert(false, ExcNotImplemented());
                            }
                        }
                    else
                      {
                        if (integrate == false)
                          for (unsigned int i = 0; i < dofs_per_face; ++i)
                            temp1[i] = VectorizedArrayType();

                        for (unsigned int v = 0; v < n_filled_lanes; ++v)
                          for (unsigned int i = 0; i < dofs_per_face; ++i)
                            proc.value(
                              temp1[reorientate(
                                n_face_orientations == 1 ? 0 : v, i)][v],
                              global_vector_ptr
                                [indices[v] +
                                 index_array_nodal
                                     [n_face_orientations == 1 ? 0 : v][i] *
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
                const unsigned int *indices =
                  &dof_info.dof_indices_contiguous[dof_access_index]
                                                  [cell *
                                                   VectorizedArrayType::size()];
                Number2_ *vector_ptr =
                  global_vector_ptr + comp * static_dofs_per_component +
                  dof_info
                    .component_dof_indices_offset[active_fe_index]
                                                 [first_selected_component];

                const unsigned int n_filled_lanes =
                  dof_info.n_vectorization_lanes_filled[dof_access_index][cell];

                const bool vectorization_possible =
                  (n_face_orientations == 1) &&
                  (n_filled_lanes == VectorizedArrayType::size()) &&
                  (sm_ptr != nullptr);

                std::array<Number2_ *, VectorizedArrayType::size()>
                  vector_ptrs = {};

                if (vectorization_possible == false)
                  {
                    if (n_face_orientations == 1)
                      {
                        for (unsigned int v = 0; v < n_filled_lanes; ++v)
                          if (sm_ptr == nullptr)
                            {
                              vector_ptrs[v] = vector_ptr + indices[v];
                            }
                          else
                            {
                              const auto &temp =
                                dof_info.dof_indices_contiguous_sm
                                  [dof_access_index]
                                  [cell * VectorizedArrayType::size() + v];
                              vector_ptrs[v] = const_cast<Number *>(
                                sm_ptr->operator[](temp.first).data() +
                                temp.second + comp * static_dofs_per_component +
                                dof_info.component_dof_indices_offset
                                  [active_fe_index][first_selected_component]);
                            }
                      }
                    else if (n_face_orientations == VectorizedArrayType::size())
                      {
                        for (unsigned int v = 0;
                             v < VectorizedArrayType::size();
                             ++v)
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
                                    dof_info.dof_indices_contiguous_sm
                                      [dof_access_index][cells[v]];
                                  vector_ptrs[v] = const_cast<Number *>(
                                    sm_ptr->operator[](temp.first).data() +
                                    temp.second +
                                    comp * static_dofs_per_component +
                                    dof_info.component_dof_indices_offset
                                      [active_fe_index]
                                      [first_selected_component]);
                                }
                            }
                      }
                    else
                      {
                        Assert(false, ExcNotImplemented());
                      }
                  }

                if (do_gradients == true &&
                    data.element_type ==
                      MatrixFreeFunctions::tensor_symmetric_hermite)
                  {
                    if (vectorization_possible)
                      for (unsigned int i = 0; i < dofs_per_face; ++i)
                        {
                          const unsigned int ind1 =
                            index_array_hermite[0][2 * i];
                          const unsigned int ind2 =
                            index_array_hermite[0][2 * i + 1];
                          const unsigned int i_ = reorientate(0, i);

                          proc.hermite_grad_vectorized_indexed(
                            temp1[i_],
                            temp1[i_ + dofs_per_face],
                            vector_ptr + ind1,
                            vector_ptr + ind2,
                            grad_weight,
                            indices,
                            indices);
                        }
                    else if (n_face_orientations == 1)
                      for (unsigned int i = 0; i < dofs_per_face; ++i)
                        {
                          const unsigned int ind1 =
                            index_array_hermite[0][2 * i];
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
                            for (unsigned int v = n_filled_lanes;
                                 v < VectorizedArrayType::size();
                                 ++v)
                              {
                                temp1[i_][v]                 = 0.0;
                                temp1[i_ + dofs_per_face][v] = 0.0;
                              }
                        }
                    else
                      {
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
                                                        indices);
                        }
                    else if (n_face_orientations == 1)
                      for (unsigned int i = 0; i < dofs_per_face; ++i)
                        {
                          const unsigned int ind = index_array_nodal[0][i];
                          const unsigned int i_  = reorientate(0, i);

                          for (unsigned int v = 0; v < n_filled_lanes; ++v)
                            proc.value(temp1[i_][v], vector_ptrs[v][ind]);

                          if (integrate == false)
                            for (unsigned int v = n_filled_lanes;
                                 v < VectorizedArrayType::size();
                                 ++v)
                              temp1[i_][v] = 0.0;
                        }
                    else
                      for (unsigned int i = 0; i < dofs_per_face; ++i)
                        {
                          for (unsigned int v = 0;
                               v < VectorizedArrayType::size();
                               ++v)
                            if (cells[v] != numbers::invalid_unsigned_int)
                              proc.value(
                                temp1[reorientate(v, i)][v],
                                vector_ptrs[v][index_array_nodal[v][i]]);
                        }
                  }
              }
            else
              {
                // case 5: default vector access
                // for the integrate_scatter path (integrate == true), we
                // need to only prepare the data in this function for all
                // components to later call distribute_local_to_global();
                // for the gather_evaluate path (integrate == false), we
                // instead want to leave early because we need to get the
                // vector data from somewhere else
                proc.default_operation(temp1, comp);
                if (integrate)
                  accesses_global_vector = false;
                else
                  return false;
              }
          }
        else
          {
            // case 5: default vector access
            proc.default_operation(temp1, comp);
            if (integrate)
              accesses_global_vector = false;
            else
              return false;
          }

        if (!integrate)
          proc.in_face_operation(temp1, comp);
      }

    if (!integrate &&
        (face_orientations[0] > 0 &&
         subface_index < GeometryInfo<dim>::max_children_per_cell))
      {
        AssertDimension(n_face_orientations, 1);
        adjust_for_face_orientation(dim,
                                    n_components,
                                    face_orientations[0],
                                    orientation_map,
                                    false,
                                    do_values,
                                    do_gradients,
                                    data.n_q_points_face,
                                    scratch_data,
                                    values_quad,
                                    gradients_quad);
      }

    return accesses_global_vector;
  }



  template <int dim,
            typename Number,
            typename VectorizedArrayType,
            typename Number2 = Number>
  struct FEFaceEvaluationImplGatherEvaluateSelector
  {
    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                          n_components,
        const unsigned int                          n_face_orientations,
        const Number2 *                             src_ptr,
        const std::vector<ArrayView<const Number>> *sm_ptr,
        const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
        const MatrixFreeFunctions::DoFInfo &                       dof_info,
        VectorizedArrayType *                                      values_quad,
        VectorizedArrayType *gradients_quad,
        VectorizedArrayType *scratch_data,
        const bool           evaluate_values,
        const bool           evaluate_gradients,
        const unsigned int   active_fe_index,
        const unsigned int   first_selected_component,
        const std::array<unsigned int, VectorizedArrayType::size()> cells,
        const std::array<unsigned int, VectorizedArrayType::size()> face_nos,
        const unsigned int                                 subface_index,
        const MatrixFreeFunctions::DoFInfo::DoFAccessIndex dof_access_index,
        const std::array<unsigned int, VectorizedArrayType::size()>
                                      face_orientations,
        const Table<2, unsigned int> &orientation_map)
    {
      if (src_ptr == nullptr)
        return false;

      if (data.element_type == MatrixFreeFunctions::tensor_none)
        return false;

      (void)sm_ptr;

      Processor<fe_degree, n_q_points_1d> p(n_components,
                                            false,
                                            src_ptr,
                                            sm_ptr,
                                            data,
                                            dof_info,
                                            values_quad,
                                            gradients_quad,
                                            scratch_data,
                                            evaluate_values,
                                            evaluate_gradients,
                                            active_fe_index,
                                            first_selected_component,
                                            cells,
                                            face_nos,
                                            subface_index,
                                            dof_access_index,
                                            face_orientations,
                                            orientation_map);

      if (n_face_orientations == VectorizedArrayType::size())
        return fe_face_evaluation_process_and_io<VectorizedArrayType::size()>(
          p);
      else
        return fe_face_evaluation_process_and_io<1>(p);
    }

  private:
    template <int fe_degree, int n_q_points_1d>
    struct Processor
    {
      static const int dim_           = dim;
      static const int fe_degree_     = fe_degree;
      static const int n_q_points_1d_ = n_q_points_1d;
      using VectorizedArrayType_      = VectorizedArrayType;
      using Number_                   = Number;
      using Number2_                  = const Number2;

      Processor(
        const unsigned int                          n_components,
        const bool                                  integrate,
        const Number2 *                             global_vector_ptr,
        const std::vector<ArrayView<const Number>> *sm_ptr,
        const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
        const MatrixFreeFunctions::DoFInfo &                       dof_info,
        VectorizedArrayType *                                      values_quad,
        VectorizedArrayType *gradients_quad,
        VectorizedArrayType *scratch_data,
        const bool           do_values,
        const bool           do_gradients,
        const unsigned int   active_fe_index,
        const unsigned int   first_selected_component,
        const std::array<unsigned int, VectorizedArrayType::size()> cells,
        const std::array<unsigned int, VectorizedArrayType::size()> face_nos,
        const unsigned int                                 subface_index,
        const MatrixFreeFunctions::DoFInfo::DoFAccessIndex dof_access_index,
        const std::array<unsigned int, VectorizedArrayType::size()>
                                      face_orientations,
        const Table<2, unsigned int> &orientation_map)
        : n_components(n_components)
        , integrate(integrate)
        , global_vector_ptr(global_vector_ptr)
        , sm_ptr(sm_ptr)
        , data(data)
        , dof_info(dof_info)
        , values_quad(values_quad)
        , gradients_quad(gradients_quad)
        , scratch_data(scratch_data)
        , do_values(do_values)
        , do_gradients(do_gradients)
        , active_fe_index(active_fe_index)
        , first_selected_component(first_selected_component)
        , cells(cells)
        , face_nos(face_nos)
        , subface_index(subface_index)
        , dof_access_index(dof_access_index)
        , face_orientations(face_orientations)
        , orientation_map(orientation_map)
      {}

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
                   const T2 &src_ptr_2,
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

      template <typename T1>
      void
      default_operation(const T1 &, const unsigned int)
      {
        // case 5)
      }

      template <typename T1>
      void
      in_face_operation(T1 &temp1, const unsigned int comp)
      {
        const unsigned int dofs_per_face =
          fe_degree > -1 ?
            Utilities::pow(fe_degree + 1, dim - 1) :
            Utilities::pow(data.data.front().fe_degree + 1, dim - 1);
        const unsigned int n_q_points =
          fe_degree > -1 ? Utilities::pow(n_q_points_1d, dim - 1) :
                           data.n_q_points_face;
        if (fe_degree > -1 &&
            subface_index >= GeometryInfo<dim>::max_children_per_cell &&
            data.element_type <= MatrixFreeFunctions::tensor_symmetric)
          FEFaceEvaluationImpl<true,
                               dim,
                               fe_degree,
                               n_q_points_1d,
                               VectorizedArrayType>::
            evaluate_in_face(/* n_components */ 1,
                             data,
                             temp1,
                             values_quad + comp * n_q_points,
                             gradients_quad + comp * dim * n_q_points,
                             scratch_data + 2 * dofs_per_face,
                             do_values,
                             do_gradients,
                             subface_index);
        else
          FEFaceEvaluationImpl<false,
                               dim,
                               fe_degree,
                               n_q_points_1d,
                               VectorizedArrayType>::
            evaluate_in_face(/* n_components */ 1,
                             data,
                             temp1,
                             values_quad + comp * n_q_points,
                             gradients_quad + comp * dim * n_q_points,
                             scratch_data + 2 * dofs_per_face,
                             do_values,
                             do_gradients,
                             subface_index);
      }

      const unsigned int                          n_components;
      const bool                                  integrate;
      const Number2 *                             global_vector_ptr;
      const std::vector<ArrayView<const Number>> *sm_ptr;
      const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data;
      const MatrixFreeFunctions::DoFInfo &                       dof_info;
      VectorizedArrayType *                                      values_quad;
      VectorizedArrayType *                                      gradients_quad;
      VectorizedArrayType *                                      scratch_data;
      const bool                                                 do_values;
      const bool                                                 do_gradients;
      const unsigned int active_fe_index;
      const unsigned int first_selected_component;
      const std::array<unsigned int, VectorizedArrayType::size()> cells;
      const std::array<unsigned int, VectorizedArrayType::size()> face_nos;
      const unsigned int                                          subface_index;
      const MatrixFreeFunctions::DoFInfo::DoFAccessIndex dof_access_index;
      const std::array<unsigned int, VectorizedArrayType::size()>
                                    face_orientations;
      const Table<2, unsigned int> &orientation_map;
    };
  };

  template <int dim,
            typename Number,
            typename VectorizedArrayType,
            typename Number2 = Number>
  struct FEFaceEvaluationImplIntegrateScatterSelector
  {
    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                           n_components,
        const unsigned int                           n_face_orientations,
        Number2 *                                    dst_ptr,
        const std::vector<ArrayView<const Number2>> *sm_ptr,
        const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
        const MatrixFreeFunctions::DoFInfo &                       dof_info,
        VectorizedArrayType *                                      values_array,
        VectorizedArrayType *                                      values_quad,
        VectorizedArrayType *gradients_quad,
        VectorizedArrayType *scratch_data,
        const bool           integrate_values,
        const bool           integrate_gradients,
        const unsigned int   active_fe_index,
        const unsigned int   first_selected_component,
        const std::array<unsigned int, VectorizedArrayType::size()> cells,
        const std::array<unsigned int, VectorizedArrayType::size()> face_nos,
        const unsigned int                                 subface_index,
        const MatrixFreeFunctions::DoFInfo::DoFAccessIndex dof_access_index,
        const std::array<unsigned int, VectorizedArrayType::size()>
                                      face_orientations,
        const Table<2, unsigned int> &orientation_map)
    {
      (void)sm_ptr;

      if (dst_ptr == nullptr ||
          data.element_type == MatrixFreeFunctions::tensor_none)
        {
          AssertDimension(n_face_orientations, 1);

          // for block vectors simply integrate
          FEFaceEvaluationImplIntegrateSelector<dim, VectorizedArrayType>::
            template run<fe_degree, n_q_points_1d>(n_components,
                                                   data,
                                                   values_array,
                                                   values_quad,
                                                   gradients_quad,
                                                   scratch_data,
                                                   integrate_values,
                                                   integrate_gradients,
                                                   face_nos[0],
                                                   subface_index,
                                                   face_orientations[0],
                                                   orientation_map);

          // default vector access
          return false;
        }


      Processor<fe_degree, n_q_points_1d> p(values_array,
                                            n_components,
                                            true,
                                            dst_ptr,
                                            sm_ptr,
                                            data,
                                            dof_info,
                                            values_quad,
                                            gradients_quad,
                                            scratch_data,
                                            integrate_values,
                                            integrate_gradients,
                                            active_fe_index,
                                            first_selected_component,
                                            cells,
                                            face_nos,
                                            subface_index,
                                            dof_access_index,
                                            face_orientations,
                                            orientation_map);

      if (n_face_orientations == VectorizedArrayType::size())
        return fe_face_evaluation_process_and_io<VectorizedArrayType::size()>(
          p);
      else
        return fe_face_evaluation_process_and_io<1>(p);
    }

  private:
    template <int fe_degree, int n_q_points_1d>
    struct Processor
    {
      static const int dim_           = dim;
      static const int fe_degree_     = fe_degree;
      static const int n_q_points_1d_ = n_q_points_1d;
      using VectorizedArrayType_      = VectorizedArrayType;
      using Number_                   = Number;
      using Number2_                  = Number2;


      Processor(
        VectorizedArrayType *                       values_array,
        const unsigned int                          n_components,
        const bool                                  integrate,
        Number2 *                                   global_vector_ptr,
        const std::vector<ArrayView<const Number>> *sm_ptr,
        const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
        const MatrixFreeFunctions::DoFInfo &                       dof_info,
        VectorizedArrayType *                                      values_quad,
        VectorizedArrayType *gradients_quad,
        VectorizedArrayType *scratch_data,
        const bool           do_values,
        const bool           do_gradients,
        const unsigned int   active_fe_index,
        const unsigned int   first_selected_component,
        const std::array<unsigned int, VectorizedArrayType::size()> cells,
        const std::array<unsigned int, VectorizedArrayType::size()> face_nos,
        const unsigned int                                 subface_index,
        const MatrixFreeFunctions::DoFInfo::DoFAccessIndex dof_access_index,
        const std::array<unsigned int, VectorizedArrayType::size()>
                                      face_orientations,
        const Table<2, unsigned int> &orientation_map)
        : values_array(values_array)
        , n_components(n_components)
        , integrate(integrate)
        , global_vector_ptr(global_vector_ptr)
        , sm_ptr(sm_ptr)
        , data(data)
        , dof_info(dof_info)
        , values_quad(values_quad)
        , gradients_quad(gradients_quad)
        , scratch_data(scratch_data)
        , do_values(do_values)
        , do_gradients(do_gradients)
        , active_fe_index(active_fe_index)
        , first_selected_component(first_selected_component)
        , cells(cells)
        , face_nos(face_nos)
        , subface_index(subface_index)
        , dof_access_index(dof_access_index)
        , face_orientations(face_orientations)
        , orientation_map(orientation_map)
      {}

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

      template <typename T0>
      void
      default_operation(const T0 &temp1, const unsigned int comp)
      {
        // case 5: default vector access, must be handled separately, just do
        // the face-normal interpolation

        FEFaceNormalEvaluationImpl<dim, fe_degree, VectorizedArrayType>::
          template interpolate<false, false>(
            /* n_components */ 1,
            data,
            temp1,
            values_array + comp * data.dofs_per_component_on_cell,
            do_gradients,
            face_nos[0]);
      }

      template <typename T0>
      void
      in_face_operation(T0 &temp1, const unsigned int comp)
      {
        const unsigned int dofs_per_face =
          fe_degree > -1 ?
            Utilities::pow(fe_degree + 1, dim - 1) :
            Utilities::pow(data.data.front().fe_degree + 1, dim - 1);
        const unsigned int n_q_points =
          fe_degree > -1 ? Utilities::pow(n_q_points_1d, dim - 1) :
                           data.n_q_points_face;
        if (fe_degree > -1 &&
            subface_index >= GeometryInfo<dim>::max_children_per_cell &&
            data.element_type <=
              internal::MatrixFreeFunctions::tensor_symmetric)
          internal::FEFaceEvaluationImpl<true,
                                         dim,
                                         fe_degree,
                                         n_q_points_1d,
                                         VectorizedArrayType>::
            integrate_in_face(/* n_components */ 1,
                              data,
                              temp1,
                              values_quad + comp * n_q_points,
                              gradients_quad + dim * comp * n_q_points,
                              scratch_data + 2 * dofs_per_face,
                              do_values,
                              do_gradients,
                              subface_index);
        else
          internal::FEFaceEvaluationImpl<false,
                                         dim,
                                         fe_degree,
                                         n_q_points_1d,
                                         VectorizedArrayType>::
            integrate_in_face(/* n_components */ 1,
                              data,
                              temp1,
                              values_quad + comp * n_q_points,
                              gradients_quad + dim * comp * n_q_points,
                              scratch_data + 2 * dofs_per_face,
                              do_values,
                              do_gradients,
                              subface_index);
      }

      VectorizedArrayType *values_array;


      const unsigned int                          n_components;
      const bool                                  integrate;
      Number2 *                                   global_vector_ptr;
      const std::vector<ArrayView<const Number>> *sm_ptr;
      const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data;
      const MatrixFreeFunctions::DoFInfo &                       dof_info;
      VectorizedArrayType *                                      values_quad;
      VectorizedArrayType *                                      gradients_quad;
      VectorizedArrayType *                                      scratch_data;
      const bool                                                 do_values;
      const bool                                                 do_gradients;
      const unsigned int active_fe_index;
      const unsigned int first_selected_component;
      const std::array<unsigned int, VectorizedArrayType::size()> cells;
      const std::array<unsigned int, VectorizedArrayType::size()> face_nos;
      const unsigned int                                          subface_index;
      const MatrixFreeFunctions::DoFInfo::DoFAccessIndex dof_access_index;
      const std::array<unsigned int, VectorizedArrayType::size()>
                                    face_orientations;
      const Table<2, unsigned int> &orientation_map;
    };
  };



  /**
   * This struct implements the action of the inverse mass matrix operation
   * using an FEEvaluationBaseData argument.
   */
  template <int dim, typename Number>
  struct CellwiseInverseMassMatrixImplBasic
  {
    template <int fe_degree, int = 0>
    static bool
    run(const unsigned int                  n_components,
        const FEEvaluationBaseData<dim,
                                   typename Number::value_type,
                                   false,
                                   Number> &fe_eval,
        const Number *                      in_array,
        Number *                            out_array,
        typename std::enable_if<fe_degree != -1>::type * = nullptr)
    {
      constexpr unsigned int dofs_per_component =
        Utilities::pow(fe_degree + 1, dim);

      Assert(dim >= 1 || dim <= 3, ExcNotImplemented());
      Assert(fe_eval.get_shape_info().element_type <=
               MatrixFreeFunctions::tensor_symmetric,
             ExcNotImplemented());

      internal::EvaluatorTensorProduct<internal::evaluate_evenodd,
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
    run(const unsigned int                  n_components,
        const FEEvaluationBaseData<dim,
                                   typename Number::value_type,
                                   false,
                                   Number> &fe_eval,
        const Number *                      in_array,
        Number *                            out_array,
        typename std::enable_if<fe_degree == -1>::type * = nullptr)
    {
      static_assert(fe_degree == -1, "Only usable for degree -1");
      const unsigned int dofs_per_component =
        fe_eval.get_shape_info().dofs_per_component_on_cell;

      Assert(dim >= 1 || dim <= 3, ExcNotImplemented());

      internal::
        EvaluatorTensorProduct<internal::evaluate_general, dim, 0, 0, Number>
          evaluator(fe_eval.get_shape_info().data.front().inverse_shape_values,
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
   * using an FEEvaluationBaseData argument.
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

      internal::EvaluatorTensorProduct<internal::evaluate_evenodd,
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
   * using an FEEvaluationBaseData argument.
   */
  template <int dim, typename Number>
  struct CellwiseInverseMassMatrixImplTransformFromQPoints
  {
    template <int fe_degree, int = 0>
    static bool
    run(const unsigned int                  n_desired_components,
        const FEEvaluationBaseData<dim,
                                   typename Number::value_type,
                                   false,
                                   Number> &fe_eval,
        const Number *                      in_array,
        Number *                            out_array,
        typename std::enable_if<fe_degree != -1>::type * = nullptr)
    {
      const AlignedVector<Number> &inverse_shape =
        fe_eval.get_shape_info().data.front().inverse_shape_values_eo;

      constexpr unsigned int dofs_per_cell = Utilities::pow(fe_degree + 1, dim);
      internal::EvaluatorTensorProduct<internal::evaluate_evenodd,
                                       dim,
                                       fe_degree + 1,
                                       fe_degree + 1,
                                       Number>
        evaluator(AlignedVector<Number>(),
                  AlignedVector<Number>(),
                  inverse_shape);

      for (unsigned int d = 0; d < n_desired_components; ++d)
        {
          const Number *in  = in_array + d * dofs_per_cell;
          Number *      out = out_array + d * dofs_per_cell;

          if (dim == 3)
            {
              evaluator.template hessians<2, false, false>(in, out);
              evaluator.template hessians<1, false, false>(out, out);
              evaluator.template hessians<0, false, false>(out, out);
            }
          if (dim == 2)
            {
              evaluator.template hessians<1, false, false>(in, out);
              evaluator.template hessians<0, false, false>(out, out);
            }
          if (dim == 1)
            evaluator.template hessians<0, false, false>(in, out);
        }
      return false;
    }

    template <int fe_degree, int = 0>
    static bool
    run(const unsigned int                  n_desired_components,
        const FEEvaluationBaseData<dim,
                                   typename Number::value_type,
                                   false,
                                   Number> &fe_eval,
        const Number *                      in_array,
        Number *                            out_array,
        typename std::enable_if<fe_degree == -1>::type * = nullptr)
    {
      static_assert(fe_degree == -1, "Only usable for degree -1");

      const AlignedVector<Number> &inverse_shape =
        fe_eval.get_shape_info().data.front().inverse_shape_values;

      const unsigned int dofs_per_component =
        fe_eval.get_shape_info().dofs_per_component_on_cell;
      const unsigned int n_q_points = fe_eval.get_shape_info().n_q_points;

      internal::
        EvaluatorTensorProduct<internal::evaluate_general, dim, 0, 0, Number>
          evaluator(inverse_shape,
                    AlignedVector<Number>(),
                    AlignedVector<Number>(),
                    fe_eval.get_shape_info().data.front().fe_degree + 1,
                    fe_eval.get_shape_info().data.front().n_q_points_1d);

      auto temp_1 = fe_eval.get_scratch_data().begin();
      auto temp_2 = temp_1 + std::max(n_q_points, dofs_per_component);

      for (unsigned int d = 0; d < n_desired_components; ++d)
        {
          const Number *in  = in_array + d * n_q_points;
          Number *      out = out_array + d * dofs_per_component;

          if (dim == 3)
            {
              evaluator.template values<2, false, false>(in, temp_1);
              evaluator.template values<1, false, false>(temp_1, temp_2);
              evaluator.template values<0, false, false>(temp_2, out);
            }
          if (dim == 2)
            {
              evaluator.template values<1, false, false>(in, temp_1);
              evaluator.template values<0, false, false>(temp_1, out);
            }
          if (dim == 1)
            evaluator.template values<0, false, false>(in, out);
        }
      return false;
    }
  };

} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
