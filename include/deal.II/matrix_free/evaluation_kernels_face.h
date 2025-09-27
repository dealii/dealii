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


#ifndef dealii_matrix_free_evaluation_kernels_face_h
#define dealii_matrix_free_evaluation_kernels_face_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/ndarray.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/matrix_free/dof_info.h>
#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/evaluation_kernels_common.h>
#include <deal.II/matrix_free/fe_evaluation_data.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
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
    using Number2 =
      typename FEEvaluationData<dim, Number, true>::shape_info_number_type;

    using Eval = EvaluatorTensorProduct<symmetric_evaluate ? evaluate_evenodd :
                                                             evaluate_general,
                                        dim - 1,
                                        fe_degree + 1,
                                        n_q_points_1d,
                                        Number,
                                        Number2>;

    static Eval
    create_evaluator_tensor_product(
      const MatrixFreeFunctions::UnivariateShapeData<Number2> &data,
      const unsigned int                                       subface_index,
      const unsigned int                                       direction)
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
      const unsigned int                                       n_components,
      const EvaluationFlags::EvaluationFlags                   evaluation_flag,
      const MatrixFreeFunctions::UnivariateShapeData<Number2> &data,
      Number                                                  *values_dofs,
      Number                                                  *values_quad,
      Number                                                  *gradients_quad,
      Number                                                  *hessians_quad,
      Number                                                  *scratch_data,
      const unsigned int                                       subface_index)
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
                  DEAL_II_NOT_IMPLEMENTED();
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
                                             Number,
                                             Number2>
                        eval_grad({}, data.shape_gradients_collocation_eo, {});
                      eval_grad.template gradients<0, true, false, 3>(
                        values_quad, gradients_quad);
                      eval_grad.template gradients<1, true, false, 3>(
                        values_quad, gradients_quad + 1);
                    }
                  else
                    {
                      // grad x
                      eval0.template gradients<0, true, false>(values_dofs,
                                                               scratch_data);
                      eval1.template values<1, true, false, 3>(scratch_data,
                                                               gradients_quad);

                      // grad y
                      eval0.template values<0, true, false>(values_dofs,
                                                            scratch_data);
                      eval1.template gradients<1, true, false, 3>(
                        scratch_data, gradients_quad + 1);

                      if ((evaluation_flag & EvaluationFlags::values) != 0u)
                        eval1.template values<1, true, false>(scratch_data,
                                                              values_quad);
                    }
                  // grad z
                  eval0.template values<0, true, false>(values_dofs + n_dofs,
                                                        scratch_data);
                  eval1.template values<1, true, false, 3>(scratch_data,
                                                           gradients_quad + 2);

                  break;
                case 2:
                  eval0.template values<0, true, false, 2>(values_dofs + n_dofs,
                                                           gradients_quad + 1);
                  eval0.template gradients<0, true, false, 2>(values_dofs,
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
      const unsigned int                                       n_components,
      const EvaluationFlags::EvaluationFlags                   integration_flag,
      const MatrixFreeFunctions::UnivariateShapeData<Number2> &data,
      Number                                                  *values_dofs,
      Number                                                  *values_quad,
      Number                                                  *gradients_quad,
      Number                                                  *hessians_quad,
      Number                                                  *scratch_data,
      const unsigned int                                       subface_index)
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
                  DEAL_II_NOT_IMPLEMENTED();
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
                  eval1.template values<1, false, false, 3>(gradients_quad + 2,
                                                            scratch_data);
                  eval0.template values<0, false, false>(scratch_data,
                                                         values_dofs + n_dofs);
                  if (symmetric_evaluate &&
                      use_collocation_evaluation(fe_degree, n_q_points_1d))
                    {
                      EvaluatorTensorProduct<evaluate_evenodd,
                                             dim - 1,
                                             n_q_points_1d,
                                             n_q_points_1d,
                                             Number,
                                             Number2>
                        eval_grad({}, data.shape_gradients_collocation_eo, {});
                      if ((integration_flag & EvaluationFlags::values) != 0u)
                        eval_grad.template gradients<1, false, true, 3>(
                          gradients_quad + 1, values_quad);
                      else
                        eval_grad.template gradients<1, false, false, 3>(
                          gradients_quad + 1, values_quad);
                      eval_grad.template gradients<0, false, true, 3>(
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
                          eval1.template gradients<1, false, true, 3>(
                            gradients_quad + 1, scratch_data);
                        }
                      else
                        eval1.template gradients<1, false, false, 3>(
                          gradients_quad + 1, scratch_data);

                      // grad y
                      eval0.template values<0, false, false>(scratch_data,
                                                             values_dofs);

                      // grad x
                      eval1.template values<1, false, false, 3>(gradients_quad,
                                                                scratch_data);
                      eval0.template gradients<0, false, true>(scratch_data,
                                                               values_dofs);
                    }
                  break;
                case 2:
                  eval0.template values<0, false, false, 2>(gradients_quad + 1,
                                                            values_dofs +
                                                              n_dofs);
                  eval0.template gradients<0, false, false, 2>(gradients_quad,
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
    using Number2 =
      typename FEEvaluationData<dim, Number, true>::shape_info_number_type;

    /**
     * Apply the sum factorization kernels within the face for Raviart-Thomas
     * elements for either evaluation or integration
     */
    template <bool do_integrate>
    static inline void
    evaluate_or_integrate_in_face(
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      const std::vector<MatrixFreeFunctions::UnivariateShapeData<Number2>>
                        &shape_data,
      Number            *values_dofs_in,
      Number            *values,
      Number            *gradients,
      Number            *scratch_data,
      const unsigned int subface_index,
      const unsigned int face_direction)
    {
      AssertDimension(shape_data.size(), 2);

      const int degree = fe_degree != -1 ? fe_degree : shape_data[0].fe_degree;
      const int n_rows_n = degree + 1;
      const int n_rows_t = degree;
      const dealii::ndarray<int, 3, 3> dofs_per_direction{
        {{{n_rows_n, n_rows_t, n_rows_t}},
         {{n_rows_t, n_rows_n, n_rows_t}},
         {{n_rows_t, n_rows_t, n_rows_n}}}};

      (void)scratch_data;
      (void)subface_index;
      // TODO: This is currently not implemented, but the test
      // matrix_vector_rt_face_03 apparently works without it -> check
      // if (subface_index < GeometryInfo<dim - 1>::max_children_per_cell)
      //  DEAL_II_NOT_IMPLEMENTED();

      using Eval = EvaluatorTensorProduct<evaluate_evenodd,
                                          dim - 1,
                                          (fe_degree > 0 ? fe_degree : 0),
                                          n_q_points_1d,
                                          Number,
                                          Number2>;

      std::array<int, dim> values_dofs_offsets = {};
      for (unsigned int comp = 0; comp < dim - 1; ++comp)
        {
          if (dim == 2)
            values_dofs_offsets[comp + 1] =
              values_dofs_offsets[comp] +
              3 * dofs_per_direction[comp][(face_direction + 1) % dim];
          else
            values_dofs_offsets[comp + 1] =
              values_dofs_offsets[comp] +
              3 * dofs_per_direction[comp][(face_direction + 1) % dim] *
                dofs_per_direction[comp][(face_direction + 2) % dim];
        }

      // Jacobians on faces are reordered to enable simple access with the
      // regular evaluators; to get the RT Piola transform right, we need to
      // pass through the values_dofs array in a permuted right order
      std::array<unsigned int, dim> components;
      for (unsigned int comp = 0; comp < dim; ++comp)
        components[comp] = (face_direction + comp + 1) % dim;

      for (const unsigned int comp : components)
        {
          Number *values_dofs = values_dofs_in + values_dofs_offsets[comp];

          std::array<int, 2> n_blocks{
            {dofs_per_direction[comp][(face_direction + 1) % dim],
             (dim > 2 ? dofs_per_direction[comp][(face_direction + 2) % dim] :
                        1)}};

          if constexpr (dim == 3)
            {
              EvaluatorTensorProduct<evaluate_evenodd,
                                     dim - 1,
                                     n_q_points_1d,
                                     n_q_points_1d,
                                     Number,
                                     Number2>
                eval_g({},
                       shape_data[0].shape_gradients_collocation_eo.data(),
                       {});
              if (!do_integrate)
                {
                  EvaluatorTensorProductAnisotropic<dim - 1,
                                                    fe_degree,
                                                    n_q_points_1d,
                                                    true>
                    eval;
                  // Evaluate in 3d
                  if (n_blocks[0] == n_rows_n)
                    {
                      eval.template normal<0>(shape_data[0],
                                              values_dofs,
                                              values);
                      eval.template tangential<1, 0>(shape_data[1],
                                                     values,
                                                     values);

                      if (evaluation_flag & EvaluationFlags::gradients)
                        {
                          eval.template normal<0>(shape_data[0],
                                                  values_dofs +
                                                    n_blocks[0] * n_blocks[1],
                                                  scratch_data);
                          eval.template tangential<1, 0, dim>(shape_data[1],
                                                              scratch_data,
                                                              gradients + 2);
                        }
                    }
                  else if (n_blocks[1] == n_rows_n)
                    {
                      eval.template normal<1>(shape_data[0],
                                              values_dofs,
                                              values);
                      eval.template tangential<0, 1>(shape_data[1],
                                                     values,
                                                     values);

                      if (evaluation_flag & EvaluationFlags::gradients)
                        {
                          eval.template normal<1>(shape_data[0],
                                                  values_dofs +
                                                    n_blocks[0] * n_blocks[1],
                                                  scratch_data);
                          eval.template tangential<0, 1, dim>(shape_data[1],
                                                              scratch_data,
                                                              gradients + 2);
                        }
                    }
                  else
                    {
                      Eval eval(shape_data[1].shape_values_eo.data(), {}, {});
                      eval.template values<0, true, false>(values_dofs, values);
                      eval.template values<1, true, false>(values, values);
                      if (evaluation_flag & EvaluationFlags::gradients)
                        {
                          eval.template values<0, true, false>(values_dofs +
                                                                 n_blocks[0] *
                                                                   n_blocks[1],
                                                               scratch_data);
                          eval.template values<1, true, false, dim>(
                            scratch_data, gradients + 2);
                        }
                    }
                  if (evaluation_flag & EvaluationFlags::gradients)
                    {
                      eval_g.template gradients<0, true, false, dim>(values,
                                                                     gradients);
                      eval_g.template gradients<1, true, false, dim>(values,
                                                                     gradients +
                                                                       1);
                    }
                }
              else
                {
                  EvaluatorTensorProductAnisotropic<dim - 1,
                                                    fe_degree,
                                                    n_q_points_1d,
                                                    false>
                    eval;
                  // Integrate in 3d
                  if (evaluation_flag & EvaluationFlags::gradients)
                    {
                      if (evaluation_flag & EvaluationFlags::values)
                        eval_g.template gradients<0, false, true, dim>(
                          gradients, values);
                      else
                        eval_g.template gradients<0, false, false, dim>(
                          gradients, values);
                      eval_g.template gradients<1, false, true, dim>(gradients +
                                                                       1,
                                                                     values);
                    }
                  if (n_blocks[0] == n_rows_n)
                    {
                      eval.template tangential<1, 0>(shape_data[1],
                                                     values,
                                                     values);
                      eval.template normal<0>(shape_data[0],
                                              values,
                                              values_dofs);

                      if (evaluation_flag & EvaluationFlags::gradients)
                        {
                          eval.template tangential<1, 0, dim>(shape_data[1],
                                                              gradients + 2,
                                                              scratch_data);
                          eval.template normal<0>(shape_data[0],
                                                  scratch_data,
                                                  values_dofs +
                                                    n_blocks[0] * n_blocks[1]);
                        }
                    }
                  else if (n_blocks[1] == n_rows_n)
                    {
                      eval.template tangential<0, 1>(shape_data[1],
                                                     values,
                                                     values);
                      eval.template normal<1>(shape_data[0],
                                              values,
                                              values_dofs);

                      if (evaluation_flag & EvaluationFlags::gradients)
                        {
                          eval.template tangential<0, 1, dim>(shape_data[1],
                                                              gradients + 2,
                                                              scratch_data);
                          eval.template normal<1>(shape_data[0],
                                                  scratch_data,
                                                  values_dofs +
                                                    n_blocks[0] * n_blocks[1]);
                        }
                    }
                  else
                    {
                      Eval eval_iso(shape_data[1].shape_values_eo.data(),
                                    {},
                                    {});
                      eval_iso.template values<1, false, false>(values, values);
                      eval_iso.template values<0, false, false>(values,
                                                                values_dofs);
                      if (evaluation_flag & EvaluationFlags::gradients)
                        {
                          eval_iso.template values<1, false, false, dim>(
                            gradients + 2, scratch_data);
                          eval_iso.template values<0, false, false>(
                            scratch_data,
                            values_dofs + n_blocks[0] * n_blocks[1]);
                        }
                    }
                }
            }
          else
            {
              using EvalN = EvaluatorTensorProduct<evaluate_evenodd,
                                                   dim - 1,
                                                   fe_degree + 1,
                                                   n_q_points_1d,
                                                   Number,
                                                   Number2>;
              if (!do_integrate)
                {
                  // Evaluate in 2d
                  if (n_blocks[0] == n_rows_n)
                    {
                      EvalN eval(shape_data[0].shape_values_eo,
                                 shape_data[0].shape_gradients_eo,
                                 {});
                      eval.template values<0, true, false>(values_dofs, values);
                      if (evaluation_flag & EvaluationFlags::gradients)
                        {
                          eval.template gradients<0, true, false, dim>(
                            values_dofs, gradients);
                          eval.template values<0, true, false, dim>(
                            values_dofs + n_rows_n, gradients + 1);
                        }
                    }
                  else
                    {
                      Eval eval(shape_data[1].shape_values_eo,
                                shape_data[1].shape_gradients_eo,
                                {});
                      eval.template values<0, true, false>(values_dofs, values);
                      if (evaluation_flag & EvaluationFlags::gradients)
                        {
                          eval.template gradients<0, true, false, dim>(
                            values_dofs, gradients);
                          eval.template values<0, true, false, dim>(
                            values_dofs + n_rows_t, gradients + 1);
                        }
                    }
                }
              else
                {
                  // Integrate in 2d
                  if (n_blocks[0] == n_rows_n)
                    {
                      EvalN eval(shape_data[0].shape_values_eo,
                                 shape_data[0].shape_gradients_eo,
                                 {});
                      if (evaluation_flag & EvaluationFlags::values)
                        eval.template values<0, false, false>(values,
                                                              values_dofs);
                      if (evaluation_flag & EvaluationFlags::gradients)
                        {
                          if (evaluation_flag & EvaluationFlags::values)
                            eval.template gradients<0, false, true, dim>(
                              gradients, values_dofs);
                          else
                            eval.template gradients<0, false, false, dim>(
                              gradients, values_dofs);
                          eval.template values<0, false, false, dim>(
                            gradients + 1, values_dofs + n_rows_n);
                        }
                    }
                  else
                    {
                      Eval eval(shape_data[1].shape_values_eo,
                                shape_data[1].shape_gradients_eo,
                                {});
                      if (evaluation_flag & EvaluationFlags::values)
                        eval.template values<0, false, false>(values,
                                                              values_dofs);
                      if (evaluation_flag & EvaluationFlags::gradients)
                        {
                          if (evaluation_flag & EvaluationFlags::values)
                            eval.template gradients<0, false, true, dim>(
                              gradients, values_dofs);
                          else
                            eval.template gradients<0, false, false, dim>(
                              gradients, values_dofs);
                          eval.template values<0, false, false, dim>(
                            gradients + 1, values_dofs + n_rows_t);
                        }
                    }
                }
            }
          values += Utilities::pow(n_q_points_1d, dim - 1);
          gradients += dim * Utilities::pow(n_q_points_1d, dim - 1);
        }
    }
  };



  template <int dim, int fe_degree, typename Number>
  struct FEFaceNormalEvaluationImpl
  {
    using Number2 =
      typename FEEvaluationData<dim, Number, true>::shape_info_number_type;

    template <bool do_evaluate, bool add_into_output>
    static void
    interpolate(const unsigned int                             n_components,
                const EvaluationFlags::EvaluationFlags         flags,
                const MatrixFreeFunctions::ShapeInfo<Number2> &shape_info,
                const Number                                  *input,
                Number                                        *output,
                const unsigned int                             face_no)
    {
      Assert(static_cast<unsigned int>(fe_degree) ==
                 shape_info.data.front().fe_degree ||
               fe_degree == -1,
             ExcInternalError());
      if (shape_info.element_type == MatrixFreeFunctions::tensor_raviart_thomas)
        interpolate_raviart_thomas<do_evaluate, add_into_output>(
          n_components, input, output, flags, face_no, shape_info);
      else
        {
          const unsigned int fe_degree_ = shape_info.data.front().fe_degree;

          interpolate_generic<do_evaluate, add_into_output>(
            n_components,
            input,
            output,
            flags,
            face_no,
            fe_degree_ + 1,
            shape_info.data.front().shape_data_on_face,
            Utilities::pow(fe_degree_ + 1, dim),
            3 * Utilities::pow(fe_degree_ + 1, dim - 1));
        }
    }

    /**
     * Interpolate the values on the cell quadrature points onto a face.
     */
    template <bool do_evaluate, bool add_into_output>
    static void
    interpolate_quadrature(
      const unsigned int                             n_components,
      const EvaluationFlags::EvaluationFlags         flags,
      const MatrixFreeFunctions::ShapeInfo<Number2> &shape_info,
      const Number                                  *input,
      Number                                        *output,
      const unsigned int                             face_no)
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
                        const Number                          *input,
                        Number                                *output,
                        const EvaluationFlags::EvaluationFlags flag,
                        const unsigned int                     face_no,
                        const unsigned int                     n_points_1d,
                        const std::array<AlignedVector<Number2>, 2> &shape_data,
                        const unsigned int dofs_per_component_on_cell,
                        const unsigned int dofs_per_component_on_face)
    {
      if (face_direction == face_no / 2)
        {
          constexpr int stride_ = Utilities::pow(fe_degree + 1, face_direction);

          const int n_rows = fe_degree != -1 ? fe_degree + 1 : n_points_1d;
          const int stride = Utilities::pow(n_rows, face_direction);
          const std::array<int, 2> n_blocks{
            {(dim > 1 ? n_rows : 1), (dim > 2 ? n_rows : 1)}};
          std::array<int, 2> steps;
          if constexpr (face_direction == 0)
            steps = {{n_rows, 0}};
          else if constexpr (face_direction == 1 && dim == 2)
            steps = {{1, 0}};
          else if constexpr (face_direction == 1)
            // in 3d, the coordinate system is zx, not xz -> switch indices
            steps = {{n_rows * n_rows, -n_rows * n_rows * n_rows + 1}};
          else if constexpr (face_direction == 2)
            steps = {{1, 0}};

          for (unsigned int c = 0; c < n_components; ++c)
            {
              if (flag & EvaluationFlags::hessians)
                interpolate_to_face<fe_degree + 1,
                                    stride_,
                                    do_evaluate,
                                    add_into_output,
                                    2>(shape_data[face_no % 2].begin(),
                                       n_blocks,
                                       steps,
                                       input,
                                       output,
                                       n_rows,
                                       stride);
              else if (flag & EvaluationFlags::gradients)
                interpolate_to_face<fe_degree + 1,
                                    stride_,
                                    do_evaluate,
                                    add_into_output,
                                    1>(shape_data[face_no % 2].begin(),
                                       n_blocks,
                                       steps,
                                       input,
                                       output,
                                       n_rows,
                                       stride);
              else
                interpolate_to_face<fe_degree + 1,
                                    stride_,
                                    do_evaluate,
                                    add_into_output,
                                    0>(shape_data[face_no % 2].begin(),
                                       n_blocks,
                                       steps,
                                       input,
                                       output,
                                       n_rows,
                                       stride);
              if (do_evaluate)
                {
                  input += dofs_per_component_on_cell;
                  output += dofs_per_component_on_face;
                }
              else
                {
                  output += dofs_per_component_on_cell;
                  input += dofs_per_component_on_face;
                }
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

    template <bool do_evaluate,
              bool add_into_output,
              int  face_direction = 0,
              int  max_derivative = 0>
    static void
    interpolate_raviart_thomas(
      const unsigned int                             n_components,
      const Number                                  *input,
      Number                                        *output,
      const EvaluationFlags::EvaluationFlags         flag,
      const unsigned int                             face_no,
      const MatrixFreeFunctions::ShapeInfo<Number2> &shape_info)
    {
      if (dim == 1)
        {
          // This should never happen since the FE_RaviartThomasNodal is not
          // defined for dim = 1. It prevents compiler warnings of infinite
          // recursion.
          DEAL_II_ASSERT_UNREACHABLE();
          return;
        }

      bool increase_max_der = false;
      if ((flag & EvaluationFlags::hessians && max_derivative < 2) ||
          (flag & EvaluationFlags::gradients && max_derivative < 1))
        increase_max_der = true;

      if (face_direction == face_no / 2 && !increase_max_der)
        {
          constexpr int stride1 = Utilities::pow(fe_degree + 1, face_direction);
          constexpr int stride0 = Utilities::pow(fe_degree, face_direction);
          constexpr int stride2 = fe_degree * (fe_degree + 1);

          const int degree =
            fe_degree != -1 ? fe_degree : shape_info.data[0].fe_degree;
          const int n_rows_n = degree + 1;
          const int n_rows_t = degree;

          std::array<int, 3> strides{{1, 1, 1}};
          if (face_direction > 0)
            {
              strides[0] =
                n_rows_n * Utilities::pow(n_rows_t, face_direction - 1);
              strides[1] = n_rows_t * (face_direction == 3 ? n_rows_n : 1);
              strides[2] = Utilities::pow(n_rows_t, face_direction);
            }
          const dealii::ndarray<int, 3, 3> dofs_per_direction{
            {{{n_rows_n, n_rows_t, n_rows_t}},
             {{n_rows_t, n_rows_n, n_rows_t}},
             {{n_rows_t, n_rows_t, n_rows_n}}}};

          std::array<int, 2> steps, n_blocks;

          if constexpr (face_direction == 0)
            steps = {{degree + (face_direction == 0), 0}};
          else if constexpr (face_direction == 1 && dim == 2)
            steps = {{1, 0}};
          else if constexpr (face_direction == 1)
            // in 3d, the coordinate system is zx, not xz -> switch indices
            steps = {
              {n_rows_n * n_rows_t, -n_rows_n * n_rows_t * n_rows_t + 1}};
          else if constexpr (face_direction == 2)
            steps = {{1, 0}};

          n_blocks[0] = dofs_per_direction[0][(face_direction + 1) % dim];
          n_blocks[1] =
            dim > 2 ? dofs_per_direction[0][(face_direction + 2) % dim] : 1;

          interpolate_to_face<
            (fe_degree != -1 ? (fe_degree + (face_direction == 0)) : 0),
            ((face_direction < 2) ? stride1 : stride2),
            do_evaluate,
            add_into_output,
            max_derivative>(shape_info.data[face_direction != 0]
                              .shape_data_on_face[face_no % 2]
                              .begin(),
                            n_blocks,
                            steps,
                            input,
                            output,
                            degree + (face_direction == 0),
                            strides[0]);

          if (do_evaluate)
            {
              input += n_rows_n * Utilities::pow(n_rows_t, dim - 1);
              output += 3 * n_blocks[0] * n_blocks[1];
            }
          else
            {
              output += n_rows_n * Utilities::pow(n_rows_t, dim - 1);
              input += 3 * n_blocks[0] * n_blocks[1];
            }

          // must only change steps only for face direction 0
          if constexpr (face_direction == 0)
            steps = {{degree, 0}};

          n_blocks[0] = dofs_per_direction[1][(face_direction + 1) % dim];
          n_blocks[1] =
            dim > 2 ? dofs_per_direction[1][(face_direction + 2) % dim] : 1;

          interpolate_to_face<
            (fe_degree != -1 ? (fe_degree + (face_direction == 1)) : 0),
            ((face_direction < 2) ? stride0 : stride2),
            do_evaluate,
            add_into_output,
            max_derivative>(shape_info.data[face_direction != 1]
                              .shape_data_on_face[face_no % 2]
                              .begin(),
                            n_blocks,
                            steps,
                            input,
                            output,
                            degree + (face_direction == 1),
                            strides[1]);

          if constexpr (dim > 2)
            {
              if (do_evaluate)
                {
                  input += n_rows_n * Utilities::pow(n_rows_t, dim - 1);
                  output += 3 * n_blocks[0] * n_blocks[1];
                }
              else
                {
                  output += n_rows_n * Utilities::pow(n_rows_t, dim - 1);
                  input += 3 * n_blocks[0] * n_blocks[1];
                }

              if constexpr (face_direction == 0)
                steps = {{degree, 0}};
              else if constexpr (face_direction == 1)
                // in 3d, the coordinate system is zx, not xz -> switch indices
                steps = {
                  {n_rows_t * n_rows_t, -n_rows_n * n_rows_t * n_rows_t + 1}};
              else if constexpr (face_direction == 2)
                steps = {{1, 0}};

              n_blocks[0] = dofs_per_direction[2][(face_direction + 1) % dim];
              n_blocks[1] = dofs_per_direction[2][(face_direction + 2) % dim];

              interpolate_to_face<
                (fe_degree != -1 ? (fe_degree + (face_direction == 2)) : 0),
                stride0,
                do_evaluate,
                add_into_output,
                max_derivative>(shape_info.data[face_direction != 2]
                                  .shape_data_on_face[face_no % 2]
                                  .begin(),
                                n_blocks,
                                steps,
                                input,
                                output,
                                degree + (face_direction == 2),
                                strides[2]);
            }
        }
      else if (face_direction == face_no / 2)
        {
          // Only increase max_derivative
          interpolate_raviart_thomas<do_evaluate,
                                     add_into_output,
                                     face_direction,
                                     std::min(max_derivative + 1, 2)>(
            n_components, input, output, flag, face_no, shape_info);
        }
      else if (face_direction < dim)
        {
          if (increase_max_der)
            {
              interpolate_raviart_thomas<do_evaluate,
                                         add_into_output,
                                         std::min(face_direction + 1, dim - 1),
                                         std::min(max_derivative + 1, 2)>(
                n_components, input, output, flag, face_no, shape_info);
            }
          else
            {
              interpolate_raviart_thomas<do_evaluate,
                                         add_into_output,
                                         std::min(face_direction + 1, dim - 1),
                                         max_derivative>(
                n_components, input, output, flag, face_no, shape_info);
            }
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
  do_vectorized_gather(const Number2       *src_ptr,
                       const unsigned int  *indices,
                       VectorizedArrayType &dst)
  {
    for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
      dst[v] = src_ptr[indices[v]];
  }



  // internal helper function for reading data; specialized version where we
  // can use a dedicated gather function
  template <typename Number, std::size_t width>
  void
  do_vectorized_gather(const Number                   *src_ptr,
                       const unsigned int             *indices,
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
                            const unsigned int       *indices,
                            Number2                  *dst_ptr)
  {
    for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
      dst_ptr[indices[v]] += src[v];
  }



  // internal helper function for reading data; specialized version where we
  // can use a dedicated gather function
  template <typename Number, std::size_t width>
  void
  do_vectorized_scatter_add(const VectorizedArray<Number, width> src,
                            const unsigned int                  *indices,
                            Number                              *dst_ptr)
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
                              Number             *tmp_values,
                              Number             *values_quad,
                              Number             *gradients_quad,
                              Number             *hessians_quad)
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
                    gradients_quad[(c * n_q_points + orientation[q]) * dim + d];
              else
                for (unsigned int q = 0; q < n_q_points; ++q)
                  tmp_values[orientation[q]] =
                    gradients_quad[(c * n_q_points + q) * dim + d];
              for (unsigned int q = 0; q < n_q_points; ++q)
                gradients_quad[(c * n_q_points + q) * dim + d] = tmp_values[q];
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
    const unsigned int                    *orientation,
    const bool                             integrate,
    const std::size_t                      n_q_points,
    Number                                *tmp_values,
    VectorizedArrayType                   *values_quad,
    VectorizedArrayType                   *gradients_quad = nullptr,
    VectorizedArrayType                   *hessians_quad  = nullptr)
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
                  tmp_values[q] =
                    gradients_quad[(c * n_q_points + orientation[q]) * dim + d]
                                  [v];
              else
                for (unsigned int q = 0; q < n_q_points; ++q)
                  tmp_values[orientation[q]] =
                    gradients_quad[(c * n_q_points + q) * dim + d][v];
              for (unsigned int q = 0; q < n_q_points; ++q)
                gradients_quad[(c * n_q_points + q) * dim + d][v] =
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
    static bool
    evaluate_tensor_none(const unsigned int                     n_components,
                         const EvaluationFlags::EvaluationFlags evaluation_flag,
                         const Number                          *values_dofs,
                         FEEvaluationData<dim, Number, true>   &fe_eval)
    {
      const auto &shape_info = fe_eval.get_shape_info();
      const auto &shape_data = shape_info.data.front();
      using Number2 =
        typename FEEvaluationData<dim, Number, true>::shape_info_number_type;

      Assert((fe_eval.get_dof_access_index() ==
                MatrixFreeFunctions::DoFInfo::dof_access_cell &&
              fe_eval.is_interior_face() == false) == false,
             ExcNotImplemented());

      const unsigned int face_no          = fe_eval.get_face_no();
      const unsigned int face_orientation = fe_eval.get_face_orientation();
      const std::size_t  n_dofs     = shape_info.dofs_per_component_on_cell;
      const std::size_t  n_q_points = shape_info.n_q_points_faces[face_no];

      if (evaluation_flag & EvaluationFlags::values)
        {
          const auto *const shape_values =
            &shape_data.shape_values_face(face_no, face_orientation, 0);

          auto *out = fe_eval.begin_values();
          auto *in  = values_dofs;

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
          auto       *out = fe_eval.begin_gradients();
          const auto *in  = values_dofs;

          const auto *const shape_gradients =
            &shape_data.shape_gradients_face(face_no, face_orientation, 0);

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

      Assert(!(evaluation_flag & EvaluationFlags::hessians),
             ExcNotImplemented());

      return true;
    }

    template <int fe_degree>
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
      static void
      project_to_face(const unsigned int                     n_components,
                      const EvaluationFlags::EvaluationFlags evaluation_flag,
                      const Number                          *values_dofs,
                      FEEvaluationData<dim, Number, true>   &fe_eval,
                      const bool                             use_vectorization,
                      Number                                *temp,
                      Number                                *scratch_data)
    {
      const auto &shape_info = fe_eval.get_shape_info();

      if (use_vectorization == false)
        {
          const auto &shape_data = shape_info.data.front();

          const unsigned int dofs_per_comp_face =
            fe_degree > -1 ?
              Utilities::pow(fe_degree + 1, dim - 1) :
              Utilities::fixed_power<dim - 1>(shape_data.fe_degree + 1);
          const unsigned int dofs_per_face = n_components * dofs_per_comp_face;

          for (unsigned int v = 0; v < Number::size(); ++v)
            {
              // the loop breaks once an invalid_unsigned_int is hit for
              // all cases except the exterior faces in the ECL loop (where
              // some faces might be at the boundaries but others not)
              if (fe_eval.get_cell_ids()[v] == numbers::invalid_unsigned_int)
                {
                  for (unsigned int i = 0; i < 3 * dofs_per_face; ++i)
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

              for (unsigned int i = 0; i < 3 * dofs_per_face; ++i)
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
    }


    template <int fe_degree, int n_q_points_1d>
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
      static void
      evaluate_in_face(const unsigned int                     n_components,
                       const EvaluationFlags::EvaluationFlags evaluation_flag,
                       FEEvaluationData<dim, Number, true>   &fe_eval,
                       Number                                *temp,
                       Number                                *scratch_data)
    {
      const auto &shape_info = fe_eval.get_shape_info();
      const auto &shape_data = shape_info.data.front();

      const unsigned int     subface_index = fe_eval.get_subface_index();
      constexpr unsigned int n_q_points_1d_actual =
        fe_degree > -1 ? n_q_points_1d : 0;

      if (shape_info.element_type == MatrixFreeFunctions::tensor_raviart_thomas)
        {
          FEFaceEvaluationImplRaviartThomas<dim,
                                            fe_degree,
                                            n_q_points_1d_actual,
                                            Number>::
            template evaluate_or_integrate_in_face<false>(
              evaluation_flag,
              fe_eval.get_shape_info().data,
              temp,
              fe_eval.begin_values(),
              fe_eval.begin_gradients(),
              scratch_data,
              subface_index,
              fe_eval.get_face_no() / 2);
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
    }

#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
    static void
    adjust_quadrature_for_face_orientation(
      const unsigned int                     n_components,
      const EvaluationFlags::EvaluationFlags evaluation_flag,
      FEEvaluationData<dim, Number, true>   &fe_eval,
      const bool                             use_vectorization,
      Number                                *temp)
    {
      const auto &shape_info = fe_eval.get_shape_info();

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
                  &shape_info.face_orientations_quad(
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
          &shape_info.face_orientations_quad(fe_eval.get_face_orientation(), 0),
          false,
          shape_info.n_q_points_face,
          temp,
          fe_eval.begin_values(),
          fe_eval.begin_gradients(),
          fe_eval.begin_hessians());
    }



    template <int fe_degree, int n_q_points_1d>
    static bool
    evaluate_tensor(const unsigned int                     n_components,
                    const EvaluationFlags::EvaluationFlags evaluation_flag,
                    const Number                          *values_dofs_actual,
                    FEEvaluationData<dim, Number, true>   &fe_eval)
    {
      const auto &shape_info = fe_eval.get_shape_info();
      const auto &shape_data = shape_info.data.front();

      const unsigned int dofs_per_comp_face =
        fe_degree > -1 ?
          Utilities::pow(fe_degree + 1, dim - 1) :
          Utilities::fixed_power<dim - 1>(shape_data.fe_degree + 1);

      // Note: we always keep storage of values, 1st and 2nd derivatives in an
      // array, so reserve space for all three here
      Number *temp1 = fe_eval.get_scratch_data().begin();
      Number *temp2 = temp1 + 3 * n_components * dofs_per_comp_face;

      const Number *values_dofs =
        (shape_data.element_type == MatrixFreeFunctions::truncated_tensor) ?
          temp2 +
            2 * (std::max(fe_eval.get_shape_info().dofs_per_component_on_cell,
                          shape_info.n_q_points)) :
          values_dofs_actual;

      if (shape_data.element_type == MatrixFreeFunctions::truncated_tensor)
        embed_truncated_into_full_tensor_product<dim, fe_degree>(
          n_components,
          const_cast<Number *>(values_dofs),
          values_dofs_actual,
          fe_eval);

      bool use_vectorization = true;
      if (fe_eval.get_dof_access_index() ==
            MatrixFreeFunctions::DoFInfo::dof_access_cell &&
          fe_eval.is_interior_face() == false) // exterior faces in the ECL loop
        for (unsigned int v = 0; v < Number::size(); ++v)
          if (fe_eval.get_cell_ids()[v] != numbers::invalid_unsigned_int &&
              fe_eval.get_face_no(v) != fe_eval.get_face_no(0))
            use_vectorization = false;

      project_to_face<fe_degree>(n_components,
                                 evaluation_flag,
                                 values_dofs,
                                 fe_eval,
                                 use_vectorization,
                                 temp1,
                                 temp2);

      evaluate_in_face<fe_degree, n_q_points_1d>(
        n_components, evaluation_flag, fe_eval, temp1, temp2);

      if (dim == 3)
        adjust_quadrature_for_face_orientation(
          n_components, evaluation_flag, fe_eval, use_vectorization, temp1);

      return false;
    }

    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                     n_components,
        const EvaluationFlags::EvaluationFlags evaluation_flag,
        const Number                          *values_dofs,
        FEEvaluationData<dim, Number, true>   &fe_eval)
    {
      const auto &shape_info = fe_eval.get_shape_info();

      if (shape_info.element_type == MatrixFreeFunctions::tensor_none)
        return evaluate_tensor_none(n_components,
                                    evaluation_flag,
                                    values_dofs,
                                    fe_eval);
      else
        return evaluate_tensor<fe_degree, n_q_points_1d>(n_components,
                                                         evaluation_flag,
                                                         values_dofs,
                                                         fe_eval);
    }
  };



  template <int dim, typename Number>
  struct FEFaceEvaluationImplProjectToFaceSelector
  {
    template <int fe_degree>
    static bool
    run(const unsigned int                     n_components,
        const EvaluationFlags::EvaluationFlags evaluation_flag,
        const Number                          *values_dofs,
        FEEvaluationData<dim, Number, true>   &fe_eval)
    {
      const auto &shape_info = fe_eval.get_shape_info();
      const auto &shape_data = shape_info.data.front();

      const unsigned int dofs_per_comp_face =
        fe_degree > -1 ?
          Utilities::pow(fe_degree + 1, dim - 1) :
          Utilities::fixed_power<dim - 1>(shape_data.fe_degree + 1);

      // Note: we always keep storage of values, 1st and 2nd derivatives in an
      // array, so reserve space for all three here
      Number *temp         = fe_eval.get_scratch_data().begin();
      Number *scratch_data = temp + 3 * n_components * dofs_per_comp_face;

      bool use_vectorization = true;
      if (fe_eval.get_dof_access_index() ==
            MatrixFreeFunctions::DoFInfo::dof_access_cell &&
          fe_eval.is_interior_face() == false) // exterior faces in the ECL loop
        for (unsigned int v = 0; v < Number::size(); ++v)
          if (fe_eval.get_cell_ids()[v] != numbers::invalid_unsigned_int &&
              fe_eval.get_face_no(v) != fe_eval.get_face_no(0))
            use_vectorization = false;

      FEFaceEvaluationImplEvaluateSelector<dim, Number>::
        template project_to_face<fe_degree>(n_components,
                                            evaluation_flag,
                                            values_dofs,
                                            fe_eval,
                                            use_vectorization,
                                            temp,
                                            scratch_data);

      return false;
    }
  };



  template <int dim, typename Number>
  struct FEFaceEvaluationImplEvaluateInFaceSelector
  {
    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                     n_components,
        const EvaluationFlags::EvaluationFlags evaluation_flag,
        FEEvaluationData<dim, Number, true>   &fe_eval)
    {
      const auto &shape_info = fe_eval.get_shape_info();
      const auto &shape_data = shape_info.data.front();

      const unsigned int dofs_per_comp_face =
        fe_degree > -1 ?
          Utilities::pow(fe_degree + 1, dim - 1) :
          Utilities::fixed_power<dim - 1>(shape_data.fe_degree + 1);

      // Note: we always keep storage of values, 1st and 2nd derivatives in an
      // array, so reserve space for all three here
      Number *temp         = fe_eval.get_scratch_data().begin();
      Number *scratch_data = temp + 3 * n_components * dofs_per_comp_face;

      FEFaceEvaluationImplEvaluateSelector<dim, Number>::
        template evaluate_in_face<fe_degree, n_q_points_1d>(
          n_components, evaluation_flag, fe_eval, temp, scratch_data);

      return false;
    }
  };



  template <int dim, typename Number>
  struct FEFaceEvaluationImplIntegrateSelector
  {
    static bool
    integrate_tensor_none(
      const unsigned int                     n_components,
      const EvaluationFlags::EvaluationFlags integration_flag,
      Number                                *values_dofs,
      FEEvaluationData<dim, Number, true>   &fe_eval,
      const bool                             sum_into_values)
    {
      const auto &shape_info = fe_eval.get_shape_info();
      const auto &shape_data = shape_info.data.front();
      using Number2 =
        typename FEEvaluationData<dim, Number, true>::shape_info_number_type;

      Assert((fe_eval.get_dof_access_index() ==
                MatrixFreeFunctions::DoFInfo::dof_access_cell &&
              fe_eval.is_interior_face() == false) == false,
             ExcNotImplemented());

      const unsigned int face_no          = fe_eval.get_face_no();
      const unsigned int face_orientation = fe_eval.get_face_orientation();
      const std::size_t  n_dofs     = shape_info.dofs_per_component_on_cell;
      const std::size_t  n_q_points = shape_info.n_q_points_faces[face_no];


      if (integration_flag & EvaluationFlags::values)
        {
          const auto *const shape_values =
            &shape_data.shape_values_face(face_no, face_orientation, 0);

          auto *in  = fe_eval.begin_values();
          auto *out = values_dofs;

          for (unsigned int c = 0; c < n_components; c += 3)
            {
              if (sum_into_values)
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
              else
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
              in += 3 * n_q_points;
              out += 3 * n_dofs;
            }
        }

      if (integration_flag & EvaluationFlags::gradients)
        {
          auto *in  = fe_eval.begin_gradients();
          auto *out = values_dofs;

          const auto *const shape_gradients =
            &shape_data.shape_gradients_face(face_no, face_orientation, 0);

          for (unsigned int c = 0; c < n_components; ++c)
            {
              if (!sum_into_values &&
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
              in += 3 * n_q_points * dim;
              out += 3 * n_dofs;
            }
        }

      Assert(!(integration_flag & EvaluationFlags::hessians),
             ExcNotImplemented());

      return true;
    }

#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
    static void
    adjust_quadrature_for_face_orientation(
      const unsigned int                     n_components,
      const EvaluationFlags::EvaluationFlags integration_flag,
      FEEvaluationData<dim, Number, true>   &fe_eval,
      const bool                             use_vectorization,
      Number                                *temp)
    {
      const auto &shape_info = fe_eval.get_shape_info();

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
    }

    template <int fe_degree, int n_q_points_1d>
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
      static void
      integrate_in_face(const unsigned int                     n_components,
                        const EvaluationFlags::EvaluationFlags integration_flag,
                        FEEvaluationData<dim, Number, true>   &fe_eval,
                        Number                                *temp,
                        Number                                *scratch_data)
    {
      const auto &shape_info = fe_eval.get_shape_info();
      const auto &shape_data = shape_info.data.front();

      const unsigned int n_q_points_1d_actual =
        fe_degree > -1 ? n_q_points_1d : 0;
      const unsigned int subface_index = fe_eval.get_subface_index();

      if (shape_info.element_type == MatrixFreeFunctions::tensor_raviart_thomas)
        {
          FEFaceEvaluationImplRaviartThomas<dim,
                                            fe_degree,
                                            n_q_points_1d_actual,
                                            Number>::
            template evaluate_or_integrate_in_face<true>(
              integration_flag,
              fe_eval.get_shape_info().data,
              temp,
              fe_eval.begin_values(),
              fe_eval.begin_gradients(),
              scratch_data,
              subface_index,
              fe_eval.get_face_no() / 2);
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
    }

    template <int fe_degree>
#ifndef DEBUG
    DEAL_II_ALWAYS_INLINE
#endif
      static void
      collect_from_face(const unsigned int                     n_components,
                        const EvaluationFlags::EvaluationFlags integration_flag,
                        Number                                *values_dofs,
                        FEEvaluationData<dim, Number, true>   &fe_eval,
                        const bool    use_vectorization,
                        const Number *temp,
                        Number       *scratch_data,
                        const bool    sum_into_values)
    {
      const auto &shape_info = fe_eval.get_shape_info();
      const auto &shape_data = shape_info.data.front();

      const unsigned int dofs_per_comp_face =
        fe_degree > -1 ?
          Utilities::pow(fe_degree + 1, dim - 1) :
          Utilities::fixed_power<dim - 1>(shape_data.fe_degree + 1);
      const unsigned int dofs_per_face = n_components * dofs_per_comp_face;

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
                                                   temp,
                                                   scratch_data,
                                                   fe_eval.get_face_no(v));

              if (sum_into_values)
                for (unsigned int i = 0; i < 3 * dofs_per_face; ++i)
                  values_dofs[i][v] += scratch_data[i][v];
              else
                for (unsigned int i = 0; i < 3 * dofs_per_face; ++i)
                  values_dofs[i][v] = scratch_data[i][v];
            }
        }
      else
        {
          if (sum_into_values)
            FEFaceNormalEvaluationImpl<dim, fe_degree, Number>::
              template interpolate<false, true>(n_components,
                                                integration_flag,
                                                shape_info,
                                                temp,
                                                values_dofs,
                                                fe_eval.get_face_no());
          else
            FEFaceNormalEvaluationImpl<dim, fe_degree, Number>::
              template interpolate<false, false>(n_components,
                                                 integration_flag,
                                                 shape_info,
                                                 temp,
                                                 values_dofs,
                                                 fe_eval.get_face_no());
        }
    }

    template <int fe_degree, int n_q_points_1d>
    static bool
    integrate_tensor(const unsigned int                     n_components,
                     const EvaluationFlags::EvaluationFlags integration_flag,
                     Number                                *values_dofs_actual,
                     FEEvaluationData<dim, Number, true>   &fe_eval,
                     const bool                             sum_into_values)
    {
      const auto &shape_info = fe_eval.get_shape_info();
      const auto &shape_data = shape_info.data.front();

      const unsigned int dofs_per_comp_face =
        fe_degree > -1 ?
          Utilities::pow(fe_degree + 1, dim - 1) :
          Utilities::fixed_power<dim - 1>(shape_data.fe_degree + 1);

      Number *temp1 = fe_eval.get_scratch_data().begin();
      Number *temp2 = temp1 + 3 * n_components * dofs_per_comp_face;

      // expand dof_values to tensor product for truncated tensor products
      Number *values_dofs =
        (shape_data.element_type == MatrixFreeFunctions::truncated_tensor) ?
          temp2 + 2 * (std::max<std::size_t>(
                        fe_eval.get_shape_info().dofs_per_component_on_cell,
                        fe_eval.get_shape_info().n_q_points)) :
          values_dofs_actual;

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

      if (dim == 3)
        adjust_quadrature_for_face_orientation(
          n_components, integration_flag, fe_eval, use_vectorization, temp1);

      integrate_in_face<fe_degree, n_q_points_1d>(
        n_components, integration_flag, fe_eval, temp1, temp2);

      collect_from_face<fe_degree>(n_components,
                                   integration_flag,
                                   values_dofs,
                                   fe_eval,
                                   use_vectorization,
                                   temp1,
                                   temp2,
                                   sum_into_values);


      if (shape_data.element_type == MatrixFreeFunctions::truncated_tensor)
        truncate_tensor_product_to_complete_degrees<dim, fe_degree>(
          n_components, values_dofs_actual, values_dofs, fe_eval);

      return false;
    }

    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                     n_components,
        const EvaluationFlags::EvaluationFlags integration_flag,
        Number                                *values_dofs,
        FEEvaluationData<dim, Number, true>   &fe_eval,
        const bool                             sum_into_values)
    {
      const auto &shape_info = fe_eval.get_shape_info();

      if (shape_info.element_type == MatrixFreeFunctions::tensor_none)
        return integrate_tensor_none(n_components,
                                     integration_flag,
                                     values_dofs,
                                     fe_eval,
                                     sum_into_values);
      else
        return integrate_tensor<fe_degree, n_q_points_1d>(n_components,
                                                          integration_flag,
                                                          values_dofs,
                                                          fe_eval,
                                                          sum_into_values);
    }
  };



  template <int dim, typename Number>
  struct FEFaceEvaluationImplCollectFromFaceSelector
  {
    template <int fe_degree>
    static bool
    run(const unsigned int                     n_components,
        const EvaluationFlags::EvaluationFlags integration_flag,
        Number                                *values_dofs,
        FEEvaluationData<dim, Number, true>   &fe_eval,
        const bool                             sum_into_values)
    {
      const auto &shape_info = fe_eval.get_shape_info();
      const auto &shape_data = shape_info.data.front();

      const unsigned int dofs_per_comp_face =
        fe_degree > -1 ?
          Utilities::pow(fe_degree + 1, dim - 1) :
          Utilities::fixed_power<dim - 1>(shape_data.fe_degree + 1);

      Number *temp         = fe_eval.get_scratch_data().begin();
      Number *scratch_data = temp + 3 * n_components * dofs_per_comp_face;

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

      FEFaceEvaluationImplIntegrateSelector<dim, Number>::
        template collect_from_face<fe_degree>(n_components,
                                              integration_flag,
                                              values_dofs,
                                              fe_eval,
                                              use_vectorization,
                                              temp,
                                              scratch_data,
                                              sum_into_values);

      return false;
    }
  };



  template <int dim, typename Number>
  struct FEFaceEvaluationImplIntegrateInFaceSelector
  {
    template <int fe_degree, int n_q_points_1d>
    static bool
    run(const unsigned int                     n_components,
        const EvaluationFlags::EvaluationFlags integration_flag,

        FEEvaluationData<dim, Number, true> &fe_eval)
    {
      const auto &shape_info = fe_eval.get_shape_info();
      const auto &shape_data = shape_info.data.front();

      const unsigned int dofs_per_comp_face =
        fe_degree > -1 ?
          Utilities::pow(fe_degree + 1, dim - 1) :
          Utilities::fixed_power<dim - 1>(shape_data.fe_degree + 1);

      Number *temp         = fe_eval.get_scratch_data().begin();
      Number *scratch_data = temp + 3 * n_components * dofs_per_comp_face;

      FEFaceEvaluationImplIntegrateSelector<dim, Number>::
        template integrate_in_face<fe_degree, n_q_points_1d>(
          n_components, integration_flag, fe_eval, temp, scratch_data);

      return false;
    }
  };



  template <int n_face_orientations,
            typename Processor,
            typename EvaluationData,
            const bool check_face_orientations = false>
  void
  fe_face_evaluation_process_and_io(
    Processor                             &proc,
    const unsigned int                     n_components,
    const EvaluationFlags::EvaluationFlags evaluation_flag,
    typename Processor::Number2_          *global_vector_ptr,
    const std::vector<ArrayView<const typename Processor::Number2_>> *sm_ptr,
    const EvaluationData                                             &fe_eval,
    typename Processor::VectorizedArrayType_                         *temp1)
  {
    constexpr int dim         = Processor::dim_;
    constexpr int fe_degree   = Processor::fe_degree_;
    using VectorizedArrayType = typename Processor::VectorizedArrayType_;
    constexpr int n_lanes     = VectorizedArrayType::size();

    using Number   = typename Processor::Number_;
    using Number2_ = typename Processor::Number2_;

    const auto        &shape_data = fe_eval.get_shape_info().data.front();
    constexpr bool     integrate  = Processor::do_integrate;
    const unsigned int face_no    = fe_eval.get_face_no();
    const auto        &dof_info   = fe_eval.get_dof_info();
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
                  shape_data
                    .shape_data_on_face[0][fe_degree + (integrate ?
                                                          (2 - (face_no % 2)) :
                                                          (1 + (face_no % 2)))];

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

            std::array<Number2_ *, n_lanes>   vector_ptrs{{nullptr}};
            std::array<unsigned int, n_lanes> reordered_indices{
              {numbers::invalid_unsigned_int}};

            if (vectorization_possible == false)
              {
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
                    DEAL_II_NOT_IMPLEMENTED();
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
                else
                  {
                    if constexpr (n_face_orientations == 1)
                      for (unsigned int i = 0; i < dofs_per_face; ++i)
                        {
                          const unsigned int ind = index_array_nodal[0][i];
                          const unsigned int i_  = reorientate(0, i);

                          for (unsigned int v = 0; v < n_filled_lanes; ++v)
                            proc.value(temp1[i_][v], vector_ptrs[v][ind]);

                          if constexpr (integrate == false)
                            for (unsigned int v = n_filled_lanes; v < n_lanes;
                                 ++v)
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
          }
        else
          {
            // We should not end up here, this should be caught by
            // FEFaceEvaluationImplGatherEvaluateSelector::supports()
            DEAL_II_ASSERT_UNREACHABLE();
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
        const Number2                                    *src_ptr,
        const std::vector<ArrayView<const Number2>>      *sm_ptr,
        FEEvaluationData<dim, VectorizedArrayType, true> &fe_eval)
    {
      Assert(fe_degree > -1, ExcInternalError());
      Assert(fe_eval.get_shape_info().element_type <=
               MatrixFreeFunctions::tensor_symmetric_no_collocation,
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

    template <typename Number3>
    static bool
    supports(const EvaluationFlags::EvaluationFlags             evaluation_flag,
             const MatrixFreeFunctions::ShapeInfo<Number3>     &shape_info,
             const Number2                                     *vector_ptr,
             MatrixFreeFunctions::DoFInfo::IndexStorageVariants storage)
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
            MatrixFreeFunctions::tensor_symmetric_no_collocation ||
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
      hermite_grad_vectorized(T0       &temp_1,
                              T0       &temp_2,
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
      hermite_grad_vectorized_indexed(T0       &temp_1,
                                      T0       &temp_2,
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
      hermite_grad(T0       &temp_1,
                   T0       &temp_2,
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
        Number2                                          *dst_ptr,
        const std::vector<ArrayView<const Number2>>      *sm_ptr,
        FEEvaluationData<dim, VectorizedArrayType, true> &fe_eval)
    {
      Assert(fe_degree > -1, ExcInternalError());
      Assert(fe_eval.get_shape_info().element_type <=
               MatrixFreeFunctions::tensor_symmetric_no_collocation,
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
                   T1       &dst_ptr_1,
                   T1       &dst_ptr_2,
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
} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
