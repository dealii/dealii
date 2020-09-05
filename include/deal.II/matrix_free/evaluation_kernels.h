// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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
    Number *temp1;
    Number *temp2;
    if (temp_size == 0)
      {
        temp1 = scratch_data;
        temp2 = temp1 + std::max(Utilities::fixed_power<dim>(
                                   shape_info.data.front().fe_degree + 1),
                                 Utilities::fixed_power<dim>(
                                   shape_info.data.front().n_q_points_1d));
      }
    else
      {
        temp1 = scratch_data;
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
    Number *temp1;
    Number *temp2;
    if (temp_size == 0)
      {
        temp1 = scratch_data;
        temp2 = temp1 + std::max(Utilities::fixed_power<dim>(
                                   shape_info.data.front().fe_degree + 1),
                                 Utilities::fixed_power<dim>(
                                   shape_info.data.front().n_q_points_1d));
      }
    else
      {
        temp1 = scratch_data;
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
     * @param transformation_matrix The coefficient matrix handed in as a
     *                     vector, using @p basis_size_1 rows and @p basis_size_2
     *                     columns if interpreted as a matrix.
     * @param coefficients The array of coefficients by which the result is
     *                     multiplied. Its length must be either
     *                     basis_size_2^dim or n_components*basis_size_2^dim
     * @param values_in    The array of the input of size basis_size_2^dim. It
     *                     may alias with values_out
     * @param scratch_data Array to hold temporary data during the operation.
     *                     Must be of length basis_size_2^dim
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



  template <bool symmetric_evaluate,
            int  dim,
            int  fe_degree,
            int  n_q_points_1d,
            int  n_components,
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
    evaluate_in_face(const MatrixFreeFunctions::ShapeInfo<Number> &data,
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
    integrate_in_face(const MatrixFreeFunctions::ShapeInfo<Number> &data,
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



  template <int dim, int fe_degree, int n_components, typename Number>
  struct FEFaceNormalEvaluationImpl
  {
    template <bool do_evaluate, bool add_into_output>
    static void
    interpolate(const MatrixFreeFunctions::ShapeInfo<Number> &data,
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
    interpolate_quadrature(const MatrixFreeFunctions::ShapeInfo<Number> &data,
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
    interpolate_generic(const Number *     input,
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



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType,
            typename Number2 = Number>
  struct FEFaceEvaluationSelector
  {
    static void
    evaluate(const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
             const VectorizedArrayType *   values_array,
             VectorizedArrayType *         values_quad,
             VectorizedArrayType *         gradients_quad,
             VectorizedArrayType *         scratch_data,
             const bool                    evaluate_values,
             const bool                    evaluate_gradients,
             const unsigned int            face_no,
             const unsigned int            subface_index,
             const unsigned int            face_orientation,
             const Table<2, unsigned int> &orientation_map)
    {
      constexpr unsigned int static_dofs_per_face =
        fe_degree > -1 ? Utilities::pow(fe_degree + 1, dim - 1) :
                         numbers::invalid_unsigned_int;
      const unsigned int dofs_per_face =
        fe_degree > -1 ?
          static_dofs_per_face :
          Utilities::pow(data.data.front().fe_degree + 1, dim - 1);

      // we allocate small amounts of data on the stack to signal the compiler
      // that this temporary data is only needed for the calculations but the
      // final results can be discarded and need not be written back to
      // memory. For large sizes or when the dofs per face is not a compile-time
      // constant, however, we want to go to the heap in the `scratch_data`
      // variable to not risk a stack overflow.
      constexpr unsigned int stack_array_size_threshold = 100;

      VectorizedArrayType
                           temp_data[static_dofs_per_face < stack_array_size_threshold ?
                    n_components * 2ull * static_dofs_per_face :
                    1];
      VectorizedArrayType *temp1;
      if (static_dofs_per_face < stack_array_size_threshold)
        temp1 = &temp_data[0];
      else
        temp1 = scratch_data;

      FEFaceNormalEvaluationImpl<dim,
                                 fe_degree,
                                 n_components,
                                 VectorizedArrayType>::
        template interpolate<true, false>(
          data, values_array, temp1, evaluate_gradients, face_no);

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
          n_components,
          VectorizedArrayType>::evaluate_in_face(data,
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
          n_components,
          VectorizedArrayType>::evaluate_in_face(data,
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
        adjust_for_face_orientation(face_orientation,
                                    orientation_map,
                                    false,
                                    evaluate_values,
                                    evaluate_gradients,
                                    data.n_q_points_face,
                                    scratch_data,
                                    values_quad,
                                    gradients_quad);
    }

    static void
    integrate(const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
              VectorizedArrayType *         values_array,
              VectorizedArrayType *         values_quad,
              VectorizedArrayType *         gradients_quad,
              VectorizedArrayType *         scratch_data,
              const bool                    integrate_values,
              const bool                    integrate_gradients,
              const unsigned int            face_no,
              const unsigned int            subface_index,
              const unsigned int            face_orientation,
              const Table<2, unsigned int> &orientation_map)
    {
      if (face_orientation)
        adjust_for_face_orientation(face_orientation,
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

      constexpr unsigned int stack_array_size_threshold = 100;

      VectorizedArrayType
                           temp_data[static_dofs_per_face < stack_array_size_threshold ?
                    n_components * 2 * static_dofs_per_face :
                    1];
      VectorizedArrayType *temp1;
      if (static_dofs_per_face < stack_array_size_threshold)
        temp1 = &temp_data[0];
      else
        temp1 = scratch_data;

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
          n_components,
          VectorizedArrayType>::integrate_in_face(data,
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
          n_components,
          VectorizedArrayType>::integrate_in_face(data,
                                                  temp1,
                                                  values_quad,
                                                  gradients_quad,
                                                  scratch_data +
                                                    2 * n_components *
                                                      dofs_per_face,
                                                  integrate_values,
                                                  integrate_gradients,
                                                  subface_index);

      FEFaceNormalEvaluationImpl<dim,
                                 fe_degree,
                                 n_components,
                                 VectorizedArrayType>::
        template interpolate<false, false>(
          data, temp1, values_array, integrate_gradients, face_no);
    }

    template <std::size_t n_face_orientations>
    static bool
    gather_evaluate(
      const Number2 *                                            src_ptr,
      const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
      const MatrixFreeFunctions::DoFInfo &                       dof_info,
      VectorizedArrayType *                                      values_quad,
      VectorizedArrayType *                                      gradients_quad,
      VectorizedArrayType *                                      scratch_data,
      const bool         evaluate_values,
      const bool         evaluate_gradients,
      const unsigned int active_fe_index,
      const unsigned int first_selected_component,
      const std::array<unsigned int, n_face_orientations> cells,
      const std::array<unsigned int, n_face_orientations> face_nos,
      const unsigned int                                  subface_index,
      const MatrixFreeFunctions::DoFInfo::DoFAccessIndex  dof_access_index,
      const std::array<unsigned int, n_face_orientations> face_orientations,
      const Table<2, unsigned int> &                      orientation_map)
    {
      if (src_ptr == nullptr)
        {
          return false;
        }

      return process_and_io( //
        false /*=evaluate*/,
        src_ptr,
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
        orientation_map,
        [](auto &      temp_1,
           auto &      temp_2,
           const auto  src_ptr_1,
           const auto  src_ptr_2,
           const auto &grad_weight) {
          // case 1a)
          do_vectorized_read(src_ptr_1, temp_1);
          do_vectorized_read(src_ptr_2, temp_2);
          temp_2 = grad_weight * (temp_1 - temp_2);
        },
        [](auto &temp, const auto src_ptr) {
          // case 1b)
          do_vectorized_read(src_ptr, temp);
        },
        [](auto &      temp_1,
           auto &      temp_2,
           const auto  src_ptr_1,
           const auto  src_ptr_2,
           const auto &grad_weight,
           const auto &indices_1,
           const auto &indices_2) {
          // case 2a)
          do_vectorized_gather(src_ptr_1, indices_1, temp_1);
          do_vectorized_gather(src_ptr_2, indices_2, temp_2);
          temp_2 = grad_weight * (temp_1 - temp_2);
        },
        [](auto &temp, const auto src_ptr, const auto &indices) {
          // case 2b)
          do_vectorized_gather(src_ptr, indices, temp);
        },
        [](auto &      temp_1,
           auto &      temp_2,
           const auto &src_ptr_1,
           const auto &src_ptr_2,
           const auto &grad_weight) {
          // case 3a)
          temp_1 = src_ptr_1;
          temp_2 = grad_weight * (temp_1 - src_ptr_2);
        },
        [](auto &temp, const auto &src_ptr) {
          // case 3b)
          temp = src_ptr;
        },
        [&](const auto &) {
          // case 5)
        },
        [&](auto &temp1, const auto &dofs_per_face) {
          if (fe_degree > -1 &&
              subface_index >= GeometryInfo<dim>::max_children_per_cell &&
              data.element_type <= MatrixFreeFunctions::tensor_symmetric)
            FEFaceEvaluationImpl<
              true,
              dim,
              fe_degree,
              n_q_points_1d,
              n_components,
              VectorizedArrayType>::evaluate_in_face(data,
                                                     temp1,
                                                     values_quad,
                                                     gradients_quad,
                                                     scratch_data +
                                                       2 * n_components *
                                                         dofs_per_face,
                                                     evaluate_values,
                                                     evaluate_gradients,
                                                     subface_index);
          else
            FEFaceEvaluationImpl<
              false,
              dim,
              fe_degree,
              n_q_points_1d,
              n_components,
              VectorizedArrayType>::evaluate_in_face(data,
                                                     temp1,
                                                     values_quad,
                                                     gradients_quad,
                                                     scratch_data +
                                                       2 * n_components *
                                                         dofs_per_face,
                                                     evaluate_values,
                                                     evaluate_gradients,
                                                     subface_index);
        });
    }

    template <std::size_t n_face_orientations>
    static bool
    integrate_scatter(
      Number2 *                                                  dst_ptr,
      const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
      const MatrixFreeFunctions::DoFInfo &                       dof_info,
      VectorizedArrayType *                                      values_array,
      VectorizedArrayType *                                      values_quad,
      VectorizedArrayType *                                      gradients_quad,
      VectorizedArrayType *                                      scratch_data,
      const bool         integrate_values,
      const bool         integrate_gradients,
      const unsigned int active_fe_index,
      const unsigned int first_selected_component,
      const std::array<unsigned int, n_face_orientations> cells,
      const std::array<unsigned int, n_face_orientations> face_nos,
      const unsigned int                                  subface_index,
      const MatrixFreeFunctions::DoFInfo::DoFAccessIndex  dof_access_index,
      const std::array<unsigned int, n_face_orientations> face_orientations,
      const Table<2, unsigned int> &                      orientation_map)
    {
      if (dst_ptr == nullptr)
        {
          AssertDimension(face_nos.size(), 1);
          AssertDimension(face_orientations.size(), 1);

          // for block vectors simply integrate
          integrate(data,
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

      return process_and_io( //
        true /*=integrate*/,
        dst_ptr,
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
        orientation_map,
        [](const auto &temp_1,
           const auto &temp_2,
           auto        dst_ptr_1,
           auto        dst_ptr_2,
           const auto &grad_weight) {
          // case 1a)
          const VectorizedArrayType val  = temp_1 - grad_weight * temp_2;
          const VectorizedArrayType grad = grad_weight * temp_2;
          do_vectorized_add(val, dst_ptr_1);
          do_vectorized_add(grad, dst_ptr_2);
        },
        [](const auto &temp, auto dst_ptr) {
          // case 1b)
          do_vectorized_add(temp, dst_ptr);
        },
        [](const auto &temp_1,
           const auto &temp_2,
           auto        dst_ptr_1,
           auto        dst_ptr_2,
           const auto &grad_weight,
           const auto &indices_1,
           const auto &indices_2) {
          // case 2a)
          const VectorizedArrayType val  = temp_1 - grad_weight * temp_2;
          const VectorizedArrayType grad = grad_weight * temp_2;
          do_vectorized_scatter_add(val, indices_1, dst_ptr_1);
          do_vectorized_scatter_add(grad, indices_2, dst_ptr_2);
        },
        [](const auto &temp, auto dst_ptr, const auto &indices) {
          // case 2b)
          do_vectorized_scatter_add(temp, indices, dst_ptr);
        },
        [](const auto &temp_1,
           const auto &temp_2,
           auto &      dst_ptr_1,
           auto &      dst_ptr_2,
           const auto &grad_weight) {
          // case 3a)
          const Number val  = temp_1 - grad_weight * temp_2;
          const Number grad = grad_weight * temp_2;
          dst_ptr_1 += val;
          dst_ptr_2 += grad;
        },
        [](const auto &temp, auto &dst_ptr) {
          // case 3b)
          dst_ptr += temp;
        },
        [&](const auto &temp1) {
          // case 5: default vector access, must be handled separately, just do
          // the face-normal interpolation

          AssertDimension(face_nos.size(), 1);

          FEFaceNormalEvaluationImpl<dim,
                                     fe_degree,
                                     n_components,
                                     VectorizedArrayType>::
            template interpolate<false, false>(
              data, temp1, values_array, integrate_gradients, face_nos[0]);
        },
        [&](auto &temp1, const auto &dofs_per_face) {
          if (fe_degree > -1 &&
              subface_index >= GeometryInfo<dim>::max_children_per_cell &&
              data.element_type <=
                internal::MatrixFreeFunctions::tensor_symmetric)
            internal::FEFaceEvaluationImpl<
              true,
              dim,
              fe_degree,
              n_q_points_1d,
              n_components,
              VectorizedArrayType>::integrate_in_face(data,
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
            internal::FEFaceEvaluationImpl<
              false,
              dim,
              fe_degree,
              n_q_points_1d,
              n_components,
              VectorizedArrayType>::integrate_in_face(data,
                                                      temp1,
                                                      values_quad,
                                                      gradients_quad,
                                                      scratch_data +
                                                        2 * n_components *
                                                          dofs_per_face,
                                                      integrate_values,
                                                      integrate_gradients,
                                                      subface_index);
        });
    }

  private:
    template <std::size_t n_face_orientations,
              typename Number2_,
              typename Function1a,
              typename Function1b,
              typename Function2a,
              typename Function2b,
              typename Function3a,
              typename Function3b,
              typename Function5,
              typename Function0>
    static bool
    process_and_io(
      const bool integrate,
      Number2_ * global_vector_ptr,
      const MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &data,
      const MatrixFreeFunctions::DoFInfo &                       dof_info,
      VectorizedArrayType *                                      values_quad,
      VectorizedArrayType *                                      gradients_quad,
      VectorizedArrayType *                                      scratch_data,
      const bool                                                 do_values,
      const bool                                                 do_gradients,
      const unsigned int active_fe_index,
      const unsigned int first_selected_component,
      const std::array<unsigned int, n_face_orientations> cells,
      const std::array<unsigned int, n_face_orientations> face_nos,
      const unsigned int                                  subface_index,
      const MatrixFreeFunctions::DoFInfo::DoFAccessIndex  dof_access_index,
      const std::array<unsigned int, n_face_orientations> face_orientations,
      const Table<2, unsigned int> &                      orientation_map,
      const Function1a &                                  function_1a,
      const Function1b &                                  function_1b,
      const Function2a &                                  function_2a,
      const Function2b &                                  function_2b,
      const Function3a &                                  function_3a,
      const Function3b &                                  function_3b,
      const Function5 &                                   function_5,
      const Function0 &                                   function_0)
    {
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
                data.data.front().nodal_at_cell_boundaries == true) ||
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
          AssertDimension(face_orientations.size(), 1);
          adjust_for_face_orientation(face_orientations[0],
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
        fe_degree > -1 ?
          static_dofs_per_face :
          Utilities::pow(data.data.front().fe_degree + 1, dim - 1);

      constexpr unsigned int stack_array_size_threshold = 100;

      VectorizedArrayType
        temp_data[static_dofs_per_face < stack_array_size_threshold ?
                    n_components * 2 * dofs_per_face :
                    1];
      VectorizedArrayType *__restrict temp1;
      if (static_dofs_per_face < stack_array_size_threshold)
        temp1 = &temp_data[0];
      else
        temp1 = scratch_data;

      if (integrate)
        function_0(temp1, dofs_per_face);

      const unsigned int dummy = 0;

      // re-orientation
      std::array<const unsigned int *, n_face_orientations> orientation;
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
      std::array<const unsigned int *, n_face_orientations> index_array_hermite;

      if (n_face_orientations == 1)
        index_array_hermite[0] =
          (data.data.front().nodal_at_cell_boundaries == true &&
           fe_degree > 1 &&
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
      std::array<const unsigned int *, n_face_orientations> index_array_nodal;

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

      if ((do_gradients == false &&
           data.data.front().nodal_at_cell_boundaries == true) ||
          (data.element_type == MatrixFreeFunctions::tensor_symmetric_hermite &&
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
              const unsigned int dof_index =
                dof_info
                  .dof_indices_contiguous[dof_access_index]
                                         [cell * VectorizedArrayType::size()] +
                dof_info
                    .component_dof_indices_offset[active_fe_index]
                                                 [first_selected_component] *
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
                          for (unsigned int comp = 0; comp < n_components;
                               ++comp)
                            function_1a(
                              temp1[i_ + 2 * comp * dofs_per_face],
                              temp1[i_ + dofs_per_face +
                                    2 * comp * dofs_per_face],
                              global_vector_ptr + dof_index +
                                (ind1 + comp * static_dofs_per_component) *
                                  VectorizedArrayType::size(),
                              global_vector_ptr + dof_index +
                                (ind2 + comp * static_dofs_per_component) *
                                  VectorizedArrayType::size(),
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
                          for (unsigned int comp = 0; comp < n_components;
                               ++comp)
                            function_1b(temp1[i_ + 2 * comp * dofs_per_face],
                                        global_vector_ptr + dof_index +
                                          (ind +
                                           comp * static_dofs_per_component) *
                                            VectorizedArrayType::size());
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
                &dof_info
                   .dof_indices_contiguous[dof_access_index]
                                          [cell * VectorizedArrayType::size()];
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
                          for (unsigned int comp = 0; comp < n_components;
                               ++comp)
                            function_2a(temp1[i_ + 2 * comp * dofs_per_face],
                                        temp1[i_ + dofs_per_face +
                                              2 * comp * dofs_per_face],
                                        global_vector_ptr + ind1 +
                                          comp * static_dofs_per_component *
                                            VectorizedArrayType::size() +
                                          dof_info.component_dof_indices_offset
                                              [active_fe_index]
                                              [first_selected_component] *
                                            VectorizedArrayType::size(),
                                        global_vector_ptr + ind2 +
                                          comp * static_dofs_per_component *
                                            VectorizedArrayType::size() +
                                          dof_info.component_dof_indices_offset
                                              [active_fe_index]
                                              [first_selected_component] *
                                            VectorizedArrayType::size(),
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

                          const unsigned int ind = index_array_nodal[0][i] *
                                                   VectorizedArrayType::size();
                          for (unsigned int comp = 0; comp < n_components;
                               ++comp)
                            function_2b(temp1[i_ + 2 * comp * dofs_per_face],
                                        global_vector_ptr + ind +
                                          comp * static_dofs_per_component *
                                            VectorizedArrayType::size() +
                                          dof_info.component_dof_indices_offset
                                              [active_fe_index]
                                              [first_selected_component] *
                                            VectorizedArrayType::size(),
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
                    [dof_access_index][cell * VectorizedArrayType::size() + v] +
                  dof_info
                      .component_dof_indices_offset[active_fe_index]
                                                   [first_selected_component] *
                    strides[v];
              const unsigned int n_filled_lanes =
                dof_info.n_vectorization_lanes_filled[dof_access_index][cell];

              if (fe_degree > 1 && do_gradients == true)
                {
                  if (n_filled_lanes == VectorizedArrayType::size())
                    for (unsigned int comp = 0; comp < n_components; ++comp)
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
                                  (comp * static_dofs_per_component +
                                   index_array_hermite[0 /*TODO*/][2 * i]) *
                                    strides[v];
                              unsigned int ind2[VectorizedArrayType::size()];
                              DEAL_II_OPENMP_SIMD_PRAGMA
                              for (unsigned int v = 0;
                                   v < VectorizedArrayType::size();
                                   ++v)
                                ind2[v] =
                                  indices[v] +
                                  (comp * static_dofs_per_component +
                                   index_array_hermite[0 /*TODO*/][2 * i + 1]) *
                                    strides[v];
                              function_2a(temp1[i_ + 2 * comp * dofs_per_face],
                                          temp1[i_ + dofs_per_face +
                                                2 * comp * dofs_per_face],
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
                        for (unsigned int i = 0;
                             i < n_components * 2 * dofs_per_face;
                             ++i)
                          temp1[i] = VectorizedArrayType();

                      for (unsigned int v = 0; v < n_filled_lanes; ++v)
                        for (unsigned int comp = 0; comp < n_components; ++comp)
                          for (unsigned int i = 0; i < dofs_per_face; ++i)
                            {
                              const unsigned int i_ =
                                reorientate(n_face_orientations == 1 ? 0 : v,
                                            i);
                              function_3a(
                                temp1[i_ + 2 * comp * dofs_per_face][v],
                                temp1[i_ + dofs_per_face +
                                      2 * comp * dofs_per_face][v],
                                global_vector_ptr
                                  [indices[v] +
                                   (comp * static_dofs_per_component +
                                    index_array_hermite
                                      [n_face_orientations == 1 ? 0 : v]
                                      [2 * i]) *
                                     strides[v]],
                                global_vector_ptr
                                  [indices[v] +
                                   (comp * static_dofs_per_component +
                                    index_array_hermite
                                      [n_face_orientations == 1 ? 0 : v]
                                      [2 * i + 1]) *
                                     strides[v]],
                                grad_weight[n_face_orientations == 1 ? 0 : v]);
                            }
                    }
                }
              else
                {
                  if (n_filled_lanes == VectorizedArrayType::size())
                    for (unsigned int comp = 0; comp < n_components; ++comp)
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
                                         (comp * static_dofs_per_component +
                                          index_array_nodal[0][i]) *
                                           strides[v];
                              const unsigned int i_ = reorientate(0, i);
                              function_2b(temp1[i_ + 2 * comp * dofs_per_face],
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
                        for (unsigned int i = 0;
                             i < n_components * dofs_per_face;
                             ++i)
                          temp1[i] = VectorizedArrayType();

                      for (unsigned int v = 0; v < n_filled_lanes; ++v)
                        for (unsigned int comp = 0; comp < n_components; ++comp)
                          for (unsigned int i = 0; i < dofs_per_face; ++i)
                            function_3b(
                              temp1[reorientate(
                                      n_face_orientations == 1 ? 0 : v, i) +
                                    2 * comp * dofs_per_face][v],
                              global_vector_ptr
                                [indices[v] +
                                 (comp * static_dofs_per_component +
                                  index_array_nodal
                                    [n_face_orientations == 1 ? 0 : v][i]) *
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
                &dof_info
                   .dof_indices_contiguous[dof_access_index]
                                          [cell * VectorizedArrayType::size()];

              if (do_gradients == true &&
                  data.element_type ==
                    MatrixFreeFunctions::tensor_symmetric_hermite)
                {
                  for (unsigned int i = 0; i < dofs_per_face; ++i)
                    {
                      if (n_face_orientations == 1 &&
                          dof_info
                              .n_vectorization_lanes_filled[dof_access_index]
                                                           [cell] ==
                            VectorizedArrayType::size())
                        {
                          const unsigned int ind1 =
                            index_array_hermite[0][2 * i];
                          const unsigned int ind2 =
                            index_array_hermite[0][2 * i + 1];
                          const unsigned int i_ = reorientate(0, i);

                          for (unsigned int comp = 0; comp < n_components;
                               ++comp)
                            function_2a(
                              temp1[i_ + 2 * comp * dofs_per_face],
                              temp1[i_ + dofs_per_face +
                                    2 * comp * dofs_per_face],
                              global_vector_ptr +
                                comp * static_dofs_per_component + ind1 +
                                dof_info.component_dof_indices_offset
                                  [active_fe_index][first_selected_component],
                              global_vector_ptr +
                                comp * static_dofs_per_component + ind2 +
                                dof_info.component_dof_indices_offset
                                  [active_fe_index][first_selected_component],
                              grad_weight,
                              indices,
                              indices);
                        }
                      else if (n_face_orientations == 1)
                        {
                          const unsigned int ind1 =
                            index_array_hermite[0][2 * i];
                          const unsigned int ind2 =
                            index_array_hermite[0][2 * i + 1];
                          const unsigned int i_ = reorientate(0, i);

                          const unsigned int n_filled_lanes =
                            dof_info
                              .n_vectorization_lanes_filled[dof_access_index]
                                                           [cell];

                          for (unsigned int v = 0; v < n_filled_lanes; ++v)
                            for (unsigned int comp = 0; comp < n_components;
                                 ++comp)
                              function_3a(
                                temp1[i_ + 2 * comp * dofs_per_face][v],
                                temp1[i_ + dofs_per_face +
                                      2 * comp * dofs_per_face][v],
                                global_vector_ptr
                                  [comp * static_dofs_per_component + ind1 +
                                   dof_info.component_dof_indices_offset
                                     [active_fe_index]
                                     [first_selected_component] +
                                   indices[v]],
                                global_vector_ptr
                                  [comp * static_dofs_per_component + ind2 +
                                   dof_info.component_dof_indices_offset
                                     [active_fe_index]
                                     [first_selected_component] +
                                   indices[v]],
                                grad_weight[v]);

                          if (integrate == false)
                            for (unsigned int v = n_filled_lanes;
                                 v < VectorizedArrayType::size();
                                 ++v)
                              for (unsigned int comp = 0; comp < n_components;
                                   ++comp)
                                {
                                  temp1[i_ + 2 * comp * dofs_per_face][v] = 0.0;
                                  temp1[i_ + dofs_per_face +
                                        2 * comp * dofs_per_face][v]      = 0.0;
                                }
                        }
                      else
                        {
                          Assert(false, ExcNotImplemented());

                          const unsigned int n_filled_lanes =
                            dof_info
                              .n_vectorization_lanes_filled[dof_access_index]
                                                           [cell];

                          for (unsigned int v = 0; v < n_filled_lanes; ++v)
                            for (unsigned int comp = 0; comp < n_components;
                                 ++comp)
                              function_3a(
                                temp1[reorientate(v, i) +
                                      2 * comp * dofs_per_face][v],
                                temp1[reorientate(v, i) + dofs_per_face +
                                      2 * comp * dofs_per_face][v],
                                global_vector_ptr
                                  [comp * static_dofs_per_component +
                                   index_array_hermite[v][2 * i] +
                                   dof_info.component_dof_indices_offset
                                     [active_fe_index]
                                     [first_selected_component] +
                                   indices[v]],
                                global_vector_ptr
                                  [comp * static_dofs_per_component +
                                   index_array_hermite[v][2 * i + 1] +
                                   dof_info.component_dof_indices_offset
                                     [active_fe_index]
                                     [first_selected_component] +
                                   indices[v]],
                                grad_weight[v]);
                        }
                    }
                }
              else
                {
                  for (unsigned int i = 0; i < dofs_per_face; ++i)
                    for (unsigned int comp = 0; comp < n_components; ++comp)
                      {
                        if (n_face_orientations == 1 &&
                            dof_info
                                .n_vectorization_lanes_filled[dof_access_index]
                                                             [cell] ==
                              VectorizedArrayType::size())
                          {
                            const unsigned int ind = index_array_nodal[0][i];
                            const unsigned int i_  = reorientate(0, i);

                            function_2b(
                              temp1[i_ + 2 * comp * dofs_per_face],
                              global_vector_ptr +
                                comp * static_dofs_per_component + ind +
                                dof_info.component_dof_indices_offset
                                  [active_fe_index][first_selected_component],
                              indices);
                          }
                        else if (n_face_orientations == 1)
                          {
                            const unsigned int ind = index_array_nodal[0][i];
                            const unsigned int i_  = reorientate(0, i);

                            const unsigned int n_filled_lanes =
                              dof_info
                                .n_vectorization_lanes_filled[dof_access_index]
                                                             [cell];

                            for (unsigned int v = 0; v < n_filled_lanes; ++v)
                              function_3b(
                                temp1[i_ + 2 * comp * dofs_per_face][v],
                                global_vector_ptr
                                  [comp * static_dofs_per_component + ind +
                                   dof_info.component_dof_indices_offset
                                     [active_fe_index]
                                     [first_selected_component] +
                                   indices[v]]);

                            if (integrate == false)
                              for (unsigned int v = n_filled_lanes;
                                   v < VectorizedArrayType::size();
                                   ++v)
                                temp1[i_ + 2 * comp * dofs_per_face][v] = 0.0;
                          }
                        else
                          {
                            for (unsigned int v = 0;
                                 v < VectorizedArrayType::size();
                                 ++v)
                              if (cells[v] != numbers::invalid_unsigned_int)
                                function_3b(
                                  temp1[reorientate(v, i) +
                                        2 * comp * dofs_per_face][v],
                                  global_vector_ptr
                                    [index_array_nodal[v][i] +
                                     comp * static_dofs_per_component +
                                     dof_info.component_dof_indices_offset
                                       [active_fe_index]
                                       [first_selected_component] +
                                     dof_info.dof_indices_contiguous
                                       [dof_access_index][cells[v]]]);
                          }
                      }
                }
            }
          else
            {
              // case 5: default vector access
              function_5(temp1);
              return false;
            }
        }
      else
        {
          // case 5: default vector access
          function_5(temp1);
          return false;
        }

      if (!integrate)
        function_0(temp1, dofs_per_face);

      if (!integrate &&
          (face_orientations[0] > 0 &&
           subface_index < GeometryInfo<dim>::max_children_per_cell))
        {
          AssertDimension(face_orientations.size(), 1);
          adjust_for_face_orientation(face_orientations[0],
                                      orientation_map,
                                      false,
                                      do_values,
                                      do_gradients,
                                      data.n_q_points_face,
                                      scratch_data,
                                      values_quad,
                                      gradients_quad);
        }

      return true;
    }

    static void
    adjust_for_face_orientation(const unsigned int            face_orientation,
                                const Table<2, unsigned int> &orientation_map,
                                const bool                    integrate,
                                const bool                    values,
                                const bool                    gradients,
                                const unsigned int            n_q_points,
                                VectorizedArrayType *         tmp_values,
                                VectorizedArrayType *         values_quad,
                                VectorizedArrayType *         gradients_quad)
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
                    tmp_values[q] = gradients_quad[(c * dim + d) * n_q_points +
                                                   orientation[q]];
                else
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    tmp_values[orientation[q]] =
                      gradients_quad[(c * dim + d) * n_q_points + q];
                for (unsigned int q = 0; q < n_q_points; ++q)
                  gradients_quad[(c * dim + d) * n_q_points + q] =
                    tmp_values[q];
              }
        }
    }
  }; // namespace internal



  /**
   * This struct implements the action of the inverse mass matrix operation
   */
  template <int dim, int fe_degree, typename Number>
  struct CellwiseInverseMassMatrixImpl
  {
    template <typename FEEvaluationType>
    static void
    apply(const unsigned int      n_components,
          const FEEvaluationType &fe_eval,
          const Number *          in_array,
          Number *                out_array)
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
    }

    static void
    apply(const unsigned int           n_desired_components,
          const AlignedVector<Number> &inverse_shape,
          const AlignedVector<Number> &inverse_coefficients,
          const Number *               in_array,
          Number *                     out_array)
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
    }

    static void
    transform_from_q_points_to_basis(const unsigned int n_desired_components,
                                     const AlignedVector<Number> &inverse_shape,
                                     const Number *               in_array,
                                     Number *                     out_array)
    {
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
    }
  };

} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
