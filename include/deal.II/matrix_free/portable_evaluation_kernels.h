// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#ifndef dealii__evaluation_kernels_h
#define dealii__evaluation_kernels_h

#include <deal.II/base/config.h>

#include <deal.II/base/utilities.h>

#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/portable_tensor_product_kernels.h>

#include <Kokkos_Core.hpp>


DEAL_II_NAMESPACE_OPEN


namespace Portable
{
  namespace internal
  {
    /**
     * Helper function to specify whether a transformation to collocation should
     * be used: It should give correct results (first condition), we need to be
     * able to initialize the fields in shape_info.templates.h from the
     * polynomials (second condition), and it should be the most efficient
     * choice in terms of operation counts (third condition).
     */
    constexpr bool
    use_collocation_evaluation(const unsigned int fe_degree,
                               const unsigned int n_q_points_1d)
    {
      // TODO: are the conditions suit for GPU parallelization?
      return (n_q_points_1d > fe_degree) && (n_q_points_1d < 200) &&
             (n_q_points_1d <= 3 * fe_degree / 2 + 1);
    }



    /**
     * This struct performs the evaluation of function values and gradients for
     * tensor-product finite elements. There are two specialized implementation
     * classes FEEvaluationImplCollocation (for Gauss-Lobatto elements where the
     * nodal points and the quadrature points coincide and the 'values'
     * operation is identity) and FEEvaluationImplTransformToCollocation (which
     * can be transformed to a collocation space and can then use the identity
     * in these spaces), which both allow for shorter code.
     */
    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    struct FEEvaluationImpl
    {
      using TeamHandle = Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type;
      using SharedView = Kokkos::View<Number *,
                                      MemorySpace::Default::kokkos_space::
                                        execution_space::scratch_memory_space,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

      DEAL_II_HOST_DEVICE static void
      evaluate(const unsigned int                            dof_handler_index,
               const unsigned int                            n_components,
               const EvaluationFlags::EvaluationFlags        evaluation_flag,
               const typename MatrixFree<dim, Number>::Data *data)
      {
        if (evaluation_flag == EvaluationFlags::nothing)
          return;

        // the evaluator does not need temporary storage since no in-place
        // operation takes place in this function
        auto scratch_for_eval =
          Kokkos::subview(data->shared_data[dof_handler_index].scratch_pad,
                          Kokkos::make_pair(0, 0));
        EvaluatorTensorProduct<EvaluatorVariant::evaluate_general,
                               dim,
                               fe_degree + 1,
                               n_q_points_1d,
                               Number>
          eval(data->team_member,
               data->precomputed_data[dof_handler_index].shape_values,
               data->precomputed_data[dof_handler_index].shape_gradients,
               data->precomputed_data[dof_handler_index].co_shape_gradients,
               scratch_for_eval);

        for (unsigned int c = 0; c < n_components; ++c)
          {
            auto u =
              Kokkos::subview(data->shared_data[dof_handler_index].values,
                              Kokkos::ALL,
                              c);
            auto grad_u =
              Kokkos::subview(data->shared_data[dof_handler_index].gradients,
                              Kokkos::ALL,
                              Kokkos::ALL,
                              c);

            if constexpr (dim == 1)
              {
                auto temp = Kokkos::subview(
                  data->shared_data[dof_handler_index].scratch_pad,
                  Kokkos::make_pair(0, n_q_points_1d));

                if (evaluation_flag & EvaluationFlags::gradients)
                  eval.template gradients<0, true, false, false>(
                    u, Kokkos::subview(grad_u, Kokkos::ALL, 0));
                if (evaluation_flag & EvaluationFlags::values)
                  {
                    eval.template values<0, true, false, false>(u, temp);
                    populate_view<false>(data->team_member,
                                         u,
                                         temp,
                                         n_q_points_1d);
                  }
              }
            else if constexpr (dim == 2)
              {
                constexpr int temp_size = (fe_degree + 1) * n_q_points_1d;
                auto          temp      = Kokkos::subview(
                  data->shared_data[dof_handler_index].scratch_pad,
                  Kokkos::make_pair(0, temp_size));

                // grad x
                if (evaluation_flag & EvaluationFlags::gradients)
                  {
                    eval.template gradients<0, true, false, false>(u, temp);
                    eval.template values<1, true, false, false>(
                      temp, Kokkos::subview(grad_u, Kokkos::ALL, 0));
                  }

                // grad y
                eval.template values<0, true, false, false>(u, temp);
                if (evaluation_flag & EvaluationFlags::gradients)
                  eval.template gradients<1, true, false, false>(
                    temp, Kokkos::subview(grad_u, Kokkos::ALL, 1));

                // val: can use values applied in x
                if (evaluation_flag & EvaluationFlags::values)
                  eval.template values<1, true, false, false>(temp, u);
              }
            else if constexpr (dim == 3)
              {
                constexpr int temp1_size = Utilities::pow(fe_degree + 1, 2) *
                                           n_q_points_1d,
                              temp2_size = Utilities::pow(n_q_points_1d, 2) *
                                           (fe_degree + 1);

                auto temp1 = Kokkos::subview(
                  data->shared_data[dof_handler_index].scratch_pad,
                  Kokkos::make_pair(0, temp1_size));
                auto temp2 = Kokkos::subview(
                  data->shared_data[dof_handler_index].scratch_pad,
                  Kokkos::make_pair(temp1_size, temp1_size + temp2_size));

                if (evaluation_flag & EvaluationFlags::gradients)
                  {
                    // grad x
                    eval.template gradients<0, true, false, false>(u, temp1);
                    eval.template values<1, true, false, false>(temp1, temp2);
                    eval.template values<2, true, false, false>(
                      temp2, Kokkos::subview(grad_u, Kokkos::ALL, 0));
                  }

                // grad y
                eval.template values<0, true, false, false>(u, temp1);
                if (evaluation_flag & EvaluationFlags::gradients)
                  {
                    eval.template gradients<1, true, false, false>(temp1,
                                                                   temp2);
                    eval.template values<2, true, false, false>(
                      temp2, Kokkos::subview(grad_u, Kokkos::ALL, 1));
                  }

                // grad z: can use the values applied in x direction stored
                // in temp1
                eval.template values<1, true, false, false>(temp1, temp2);
                if (evaluation_flag & EvaluationFlags::gradients)
                  eval.template gradients<2, true, false, false>(
                    temp2, Kokkos::subview(grad_u, Kokkos::ALL, 2));

                // val: can use the values applied in x & y direction stored
                // in temp2
                if (evaluation_flag & EvaluationFlags::values)
                  eval.template values<2, true, false, false>(temp2, u);
              }
            else
              Assert(false, ExcMessage("dim must not exceed 3!"));
          }
      }



      DEAL_II_HOST_DEVICE static void
      integrate(const unsigned int                            dof_handler_index,
                const unsigned int                            n_components,
                const EvaluationFlags::EvaluationFlags        integration_flag,
                const typename MatrixFree<dim, Number>::Data *data)
      {
        if (integration_flag == EvaluationFlags::nothing)
          return;

        // the evaluator does not need temporary storage since no in-place
        // operation takes place in this function
        auto scratch_for_eval =
          Kokkos::subview(data->shared_data[dof_handler_index].scratch_pad,
                          Kokkos::make_pair(0, 0));
        EvaluatorTensorProduct<EvaluatorVariant::evaluate_general,
                               dim,
                               fe_degree + 1,
                               n_q_points_1d,
                               Number>
          eval(data->team_member,
               data->precomputed_data[dof_handler_index].shape_values,
               data->precomputed_data[dof_handler_index].shape_gradients,
               data->precomputed_data[dof_handler_index].co_shape_gradients,
               scratch_for_eval);

        for (unsigned int c = 0; c < n_components; ++c)
          {
            auto u =
              Kokkos::subview(data->shared_data[dof_handler_index].values,
                              Kokkos::ALL,
                              c);
            auto grad_u =
              Kokkos::subview(data->shared_data[dof_handler_index].gradients,
                              Kokkos::ALL,
                              Kokkos::ALL,
                              c);

            if constexpr (dim == 1)
              {
                auto temp = Kokkos::subview(
                  data->shared_data[dof_handler_index].scratch_pad,
                  Kokkos::make_pair(0, fe_degree + 1));

                if ((integration_flag & EvaluationFlags::values) &&
                    !(integration_flag & EvaluationFlags::gradients))
                  {
                    eval.template values<0, false, false, false>(u, temp);
                    populate_view<false>(data->team_member,
                                         u,
                                         temp,
                                         fe_degree + 1);
                  }
                if (integration_flag & EvaluationFlags::gradients)
                  {
                    if (integration_flag & EvaluationFlags::values)
                      {
                        eval.template values<0, false, false, false>(u, temp);
                        eval.template gradients<0, false, true, false>(
                          Kokkos::subview(grad_u, Kokkos::ALL, 0), temp);
                        populate_view<false>(data->team_member,
                                             u,
                                             temp,
                                             fe_degree + 1);
                      }
                    else
                      eval.template gradients<0, false, false, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
                  }
              }
            else if constexpr (dim == 2)
              {
                constexpr int temp_size = (fe_degree + 1) * n_q_points_1d;
                auto          temp      = Kokkos::subview(
                  data->shared_data[dof_handler_index].scratch_pad,
                  Kokkos::make_pair(0, temp_size));

                if ((integration_flag & EvaluationFlags::values) &&
                    !(integration_flag & EvaluationFlags::gradients))
                  {
                    eval.template values<1, false, false, false>(u, temp);
                    eval.template values<0, false, false, false>(temp, u);
                  }
                if (integration_flag & EvaluationFlags::gradients)
                  {
                    eval.template gradients<1, false, false, false>(
                      Kokkos::subview(grad_u, Kokkos::ALL, 1), temp);
                    if (integration_flag & EvaluationFlags::values)
                      eval.template values<1, false, true, false>(u, temp);
                    eval.template values<0, false, false, false>(temp, u);
                    eval.template values<1, false, false, false>(
                      Kokkos::subview(grad_u, Kokkos::ALL, 0), temp);
                    eval.template gradients<0, false, true, false>(temp, u);
                  }
              }
            else if constexpr (dim == 3)
              {
                constexpr int temp1_size = Utilities::pow(n_q_points_1d, 2) *
                                           (fe_degree + 1),
                              temp2_size = Utilities::pow(fe_degree + 1, 2) *
                                           n_q_points_1d;

                auto temp1 = Kokkos::subview(
                  data->shared_data[dof_handler_index].scratch_pad,
                  Kokkos::make_pair(0, temp1_size));
                auto temp2 = Kokkos::subview(
                  data->shared_data[dof_handler_index].scratch_pad,
                  Kokkos::make_pair(temp1_size, temp1_size + temp2_size));

                if ((integration_flag & EvaluationFlags::values) &&
                    !(integration_flag & EvaluationFlags::gradients))
                  {
                    eval.template values<2, false, false, false>(u, temp1);
                    eval.template values<1, false, false, false>(temp1, temp2);
                    eval.template values<0, false, false, false>(temp2, u);
                  }
                if (integration_flag & EvaluationFlags::gradients)
                  {
                    eval.template gradients<2, false, false, false>(
                      Kokkos::subview(grad_u, Kokkos::ALL, 2), temp1);
                    if (integration_flag & EvaluationFlags::values)
                      eval.template values<2, false, true, false>(u, temp1);
                    eval.template values<1, false, false, false>(temp1, temp2);
                    eval.template values<2, false, false, false>(
                      Kokkos::subview(grad_u, Kokkos::ALL, 1), temp1);
                    eval.template gradients<1, false, true, false>(temp1,
                                                                   temp2);
                    eval.template values<0, false, false, false>(temp2, u);
                    eval.template values<2, false, false, false>(
                      Kokkos::subview(grad_u, Kokkos::ALL, 0), temp1);
                    eval.template values<1, false, false, false>(temp1, temp2);
                    eval.template gradients<0, false, true, false>(temp2, u);
                  }
              }
            else
              Assert(false, ExcMessage("dim must not exceed 3!"));
          }
      }
    };



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
      DEAL_II_HOST_DEVICE static void
      evaluate(const unsigned int                            dof_handler_index,
               const unsigned int                            n_components,
               const EvaluationFlags::EvaluationFlags        evaluation_flag,
               const typename MatrixFree<dim, Number>::Data *data)
      {
        // since the dof values have already been stored in
        // shared_data->values, there is nothing to do if the gradients are
        // not required
        if (!(evaluation_flag & EvaluationFlags::gradients))
          return;

        constexpr int n_points = Utilities::pow(fe_degree + 1, dim);
        auto          scratch_for_eval =
          Kokkos::subview(data->shared_data[dof_handler_index].scratch_pad,
                          Kokkos::make_pair(0, n_points));

        EvaluatorTensorProduct<EvaluatorVariant::evaluate_general,
                               dim,
                               fe_degree + 1,
                               fe_degree + 1,
                               Number>
          eval(data->team_member,
               data->precomputed_data[dof_handler_index].shape_values,
               data->precomputed_data[dof_handler_index].shape_gradients,
               data->precomputed_data[dof_handler_index].co_shape_gradients,
               scratch_for_eval);

        for (unsigned int c = 0; c < n_components; ++c)
          {
            auto u =
              Kokkos::subview(data->shared_data[dof_handler_index].values,
                              Kokkos::ALL,
                              c);
            auto grad_u =
              Kokkos::subview(data->shared_data[dof_handler_index].gradients,
                              Kokkos::ALL,
                              Kokkos::ALL,
                              c);

            eval.template co_gradients<0, true, false, false>(
              u, Kokkos::subview(grad_u, Kokkos::ALL, 0));
            if constexpr (dim > 1)
              eval.template co_gradients<1, true, false, false>(
                u, Kokkos::subview(grad_u, Kokkos::ALL, 1));
            if constexpr (dim > 2)
              eval.template co_gradients<2, true, false, false>(
                u, Kokkos::subview(grad_u, Kokkos::ALL, 2));
          }
      }


      DEAL_II_HOST_DEVICE static void
      integrate(const unsigned int                            dof_handler_index,
                const unsigned int                            n_components,
                const EvaluationFlags::EvaluationFlags        integration_flag,
                const typename MatrixFree<dim, Number>::Data *data)
      {
        // since the quad values have already been stored in
        // shared_data->values, there is nothing to do if the gradients are
        // not required
        if (!(integration_flag & EvaluationFlags::gradients))
          return;

        constexpr int n_points = Utilities::pow(fe_degree + 1, dim);
        auto          scratch_for_eval =
          Kokkos::subview(data->shared_data[dof_handler_index].scratch_pad,
                          Kokkos::make_pair(0, n_points));

        EvaluatorTensorProduct<EvaluatorVariant::evaluate_general,
                               dim,
                               fe_degree + 1,
                               fe_degree + 1,
                               Number>
          eval(data->team_member,
               data->precomputed_data[dof_handler_index].shape_values,
               data->precomputed_data[dof_handler_index].shape_gradients,
               data->precomputed_data[dof_handler_index].co_shape_gradients,
               scratch_for_eval);

        for (unsigned int c = 0; c < n_components; ++c)
          {
            auto u =
              Kokkos::subview(data->shared_data[dof_handler_index].values,
                              Kokkos::ALL,
                              c);
            auto grad_u =
              Kokkos::subview(data->shared_data[dof_handler_index].gradients,
                              Kokkos::ALL,
                              Kokkos::ALL,
                              c);

            if constexpr (dim == 1)
              {
                if (integration_flag & EvaluationFlags::values)
                  eval.template co_gradients<0, false, true, false>(
                    Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
                else
                  eval.template co_gradients<0, false, false, false>(
                    Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
              }
            else if constexpr (dim == 2)
              {
                if (integration_flag & EvaluationFlags::values)
                  eval.template co_gradients<1, false, true, false>(
                    Kokkos::subview(grad_u, Kokkos::ALL, 1), u);
                else
                  eval.template co_gradients<1, false, false, false>(
                    Kokkos::subview(grad_u, Kokkos::ALL, 1), u);
                eval.template co_gradients<0, false, true, false>(
                  Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
              }
            else if constexpr (dim == 3)
              {
                if (integration_flag & EvaluationFlags::values)
                  eval.template co_gradients<2, false, true, false>(
                    Kokkos::subview(grad_u, Kokkos::ALL, 2), u);
                else
                  eval.template co_gradients<2, false, false, false>(
                    Kokkos::subview(grad_u, Kokkos::ALL, 2), u);
                eval.template co_gradients<1, false, true, false>(
                  Kokkos::subview(grad_u, Kokkos::ALL, 1), u);
                eval.template co_gradients<0, false, true, false>(
                  Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
              }
            else
              Assert(false, ExcMessage("dim must not exceed 3!"));
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
      DEAL_II_HOST_DEVICE static void
      evaluate(const unsigned int                            dof_handler_index,
               const unsigned int                            n_components,
               const EvaluationFlags::EvaluationFlags        evaluation_flag,
               const typename MatrixFree<dim, Number>::Data *data)
      {
        constexpr int scratch_size = Utilities::pow(n_q_points_1d, dim);
        auto          scratch_for_eval =
          Kokkos::subview(data->shared_data[dof_handler_index].scratch_pad,
                          Kokkos::make_pair(0, scratch_size));

        EvaluatorTensorProduct<EvaluatorVariant::evaluate_general,
                               dim,
                               fe_degree + 1,
                               n_q_points_1d,
                               Number>
          eval(data->team_member,
               data->precomputed_data[dof_handler_index].shape_values,
               data->precomputed_data[dof_handler_index].shape_gradients,
               data->precomputed_data[dof_handler_index].co_shape_gradients,
               scratch_for_eval);

        for (unsigned int c = 0; c < n_components; ++c)
          {
            auto u =
              Kokkos::subview(data->shared_data[dof_handler_index].values,
                              Kokkos::ALL,
                              c);
            auto grad_u =
              Kokkos::subview(data->shared_data[dof_handler_index].gradients,
                              Kokkos::ALL,
                              Kokkos::ALL,
                              c);

            eval.template values<0, true, false, true>(u, u);
            if constexpr (dim > 1)
              eval.template values<1, true, false, true>(u, u);
            if constexpr (dim > 2)
              eval.template values<2, true, false, true>(u, u);

            if (evaluation_flag & EvaluationFlags::gradients)
              {
                eval.template co_gradients<0, true, false, false>(
                  u, Kokkos::subview(grad_u, Kokkos::ALL, 0));
                if constexpr (dim > 1)
                  eval.template co_gradients<1, true, false, false>(
                    u, Kokkos::subview(grad_u, Kokkos::ALL, 1));
                if constexpr (dim > 2)
                  eval.template co_gradients<2, true, false, false>(
                    u, Kokkos::subview(grad_u, Kokkos::ALL, 2));
              }
          }
      }


      DEAL_II_HOST_DEVICE static void
      integrate(const unsigned int                            dof_handler_index,
                const unsigned int                            n_components,
                const EvaluationFlags::EvaluationFlags        integration_flag,
                const typename MatrixFree<dim, Number>::Data *data)
      {
        constexpr int scratch_size = Utilities::pow(n_q_points_1d, dim);
        auto          scratch_for_eval =
          Kokkos::subview(data->shared_data[dof_handler_index].scratch_pad,
                          Kokkos::make_pair(0, scratch_size));

        EvaluatorTensorProduct<EvaluatorVariant::evaluate_general,
                               dim,
                               fe_degree + 1,
                               n_q_points_1d,
                               Number>
          eval(data->team_member,
               data->precomputed_data[dof_handler_index].shape_values,
               data->precomputed_data[dof_handler_index].shape_gradients,
               data->precomputed_data[dof_handler_index].co_shape_gradients,
               scratch_for_eval);

        for (unsigned int c = 0; c < n_components; ++c)
          {
            auto u =
              Kokkos::subview(data->shared_data[dof_handler_index].values,
                              Kokkos::ALL,
                              c);
            auto grad_u =
              Kokkos::subview(data->shared_data[dof_handler_index].gradients,
                              Kokkos::ALL,
                              Kokkos::ALL,
                              c);

            // apply derivatives in collocation space
            if (integration_flag & EvaluationFlags::gradients)
              {
                if constexpr (dim == 1)
                  {
                    if (integration_flag & EvaluationFlags::values)
                      eval.template co_gradients<0, false, true, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
                    else
                      eval.template co_gradients<0, false, false, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
                  }
                else if constexpr (dim == 2)
                  {
                    if (integration_flag & EvaluationFlags::values)
                      eval.template co_gradients<1, false, true, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 1), u);
                    else
                      eval.template co_gradients<1, false, false, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 1), u);
                    eval.template co_gradients<0, false, true, false>(
                      Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
                  }
                else if constexpr (dim == 3)
                  {
                    if (integration_flag & EvaluationFlags::values)
                      eval.template co_gradients<2, false, true, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 2), u);
                    else
                      eval.template co_gradients<2, false, false, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 2), u);
                    eval.template co_gradients<1, false, true, false>(
                      Kokkos::subview(grad_u, Kokkos::ALL, 1), u);
                    eval.template co_gradients<0, false, true, false>(
                      Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
                  }
                else
                  Assert(false, ExcMessage("dim must not exceed 3!"));
              }

            // transform back to the original space
            if constexpr (dim > 2)
              eval.template values<2, false, false, true>(u, u);
            if constexpr (dim > 1)
              eval.template values<1, false, false, true>(u, u);
            eval.template values<0, false, false, true>(u, u);
          }
      }
    };

    namespace batched
    {


      /**
       * Same as Portable::MatrixFree::SharedData, but additionally contains
       * Kokkos Views for shape data.
       */
      template <int dim, typename Number>
      struct SharedData
      {
        using TeamHandle = Kokkos::TeamPolicy<
          MemorySpace::Default::kokkos_space::execution_space>::member_type;

        using SharedViewValues =
          Kokkos::View<Number **,
                       MemorySpace::Default::kokkos_space::execution_space::
                         scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
        using SharedViewGradients =
          Kokkos::View<Number ***,
                       MemorySpace::Default::kokkos_space::execution_space::
                         scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
        using SharedViewScratchPad =
          Kokkos::View<Number *,
                       MemorySpace::Default::kokkos_space::execution_space::
                         scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
        using SharedViewShapeData =
          Kokkos::View<Number *,
                       MemorySpace::Default::kokkos_space::execution_space::
                         scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

        /**
         * Allocate Kokkos Views of the right size and copy shape data from the
         * precomputed data to shared memory.
         */
        DEAL_II_HOST_DEVICE
        void
        reinit(
          const TeamHandle                                        &team_handle,
          const unsigned int                                       n_q_points,
          const typename MatrixFree<dim, Number>::PrecomputedData &gpu_data)
        {
          values      = SharedViewValues(team_handle.team_shmem(),
                                    n_q_points,
                                    gpu_data.n_components);
          gradients   = SharedViewGradients(team_handle.team_shmem(),
                                          n_q_points,
                                          dim,
                                          gpu_data.n_components);
          scratch_pad = SharedViewScratchPad(team_handle.team_shmem(),
                                             gpu_data.scratch_pad_size);

          // For symmetric collocation elements, the shape_values are the
          // identity and do not need to be stored in shared memory
          if (gpu_data.element_type >= ::dealii::internal::MatrixFreeFunctions::
                                         ElementType::tensor_symmetric)
            {
              shape_values = SharedViewShapeData(team_handle.team_shmem(),
                                                 gpu_data.shape_values.size());

              Kokkos::parallel_for(Kokkos::TeamVectorRange(
                                     team_handle, gpu_data.shape_values.size()),
                                   [&](const int i) {
                                     shape_values(i) = gpu_data.shape_values[i];
                                   });
            }

          // For symmetric elements, we need to store the shape gradients in the
          // collocation space for faster evaluation. For other elements, we
          // store the shape gradients in the reference space.
          if (gpu_data.element_type <= ::dealii::internal::MatrixFreeFunctions::
                                         ElementType::tensor_symmetric)
            {
              co_shape_gradients =
                SharedViewShapeData(team_handle.team_shmem(),
                                    gpu_data.co_shape_gradients.size());

              Kokkos::parallel_for(
                Kokkos::TeamVectorRange(team_handle,
                                        gpu_data.co_shape_gradients.size()),
                [&](const int i) {
                  co_shape_gradients(i) = gpu_data.co_shape_gradients[i];
                });
            }
          else
            {
              shape_gradients =
                SharedViewShapeData(team_handle.team_shmem(),
                                    gpu_data.shape_gradients.size());

              Kokkos::parallel_for(
                Kokkos::TeamVectorRange(team_handle,
                                        gpu_data.shape_gradients.size()),
                [&](const int i) {
                  shape_gradients(i) = gpu_data.shape_gradients[i];
                });
            }
        }

        /**
         * Memory for dof and quad values.
         */
        SharedViewValues values;

        /**
         * Memory for computed gradients in reference coordinate system.
         */
        SharedViewGradients gradients;

        /**
         * Memory for temporary arrays required by evaluation and integration.
         */
        SharedViewScratchPad scratch_pad;

        /**
         * Memory for shape values.
         */
        SharedViewShapeData shape_values;
        /**
         * Memory for shape gradients.
         */
        SharedViewShapeData shape_gradients;
        /**
         * Memory for shape gradients in the collocation space.
         */
        SharedViewShapeData co_shape_gradients;
      };


      /**
       * Same as Portable::MatrixFree::Data, but contains the above SharedData
       * struct instead of the original one. Additionally, it constains
       * information needed for batched cell loops and evaluation.
       */
      template <int dim, typename Number>
      struct BatchedData
      {
        using TeamHandle = Kokkos::TeamPolicy<
          MemorySpace::Default::kokkos_space::execution_space>::member_type;

        /**
         * TeamPolicy handle.
         */
        TeamHandle team_member;

        const unsigned int n_q_points;
        const unsigned int n_dof_handler;
        /**
         * Number of cells per batch. It is computed based on the amount of
         * shared memory requested by the kernel and the amount of shared memory
         * available on the device to saturate the device. It is not guaranteed
         * that all batches will have the same number of cells, since the last
         * batch may have fewer cells.
         */
        const int n_cells_per_batch;
        const int n_cells_in_current_batch;

        // The index of the current batch. It is retrived from the main parallel
        // loop over the batches and is used to compute the global cell index of
        // the current batch.
        const int batch_index;

        const Kokkos::Array<typename MatrixFree<dim, Number>::PrecomputedData,
                            n_max_dof_handlers> &precomputed_data;
        Kokkos::Array<SharedData<dim, Number>, n_max_dof_handlers> &shared_data;


        /**
         * Return the quadrature point index of the given cell and @p q_point index.
         * The index returned is only unique for a given MPI process.
         */
        DEAL_II_HOST_DEVICE unsigned int
        local_q_point_id(const unsigned int cell,
                         const unsigned int q_point) const
        {
          AssertIndexRange(cell, precomputed_data[0].n_cells);
          AssertIndexRange(q_point, n_q_points);

          Assert(precomputed_data[0].n_cells ==
                   precomputed_data[0].q_points.extent(1),
                 ExcInternalError("q_points array has wrong size"));
          Assert(n_q_points == precomputed_data[0].q_points.extent(0),
                 ExcInternalError("q_points array has wrong size"));

          return (precomputed_data[0].row_start /
                    precomputed_data[0].padding_length +
                  cell) *
                   n_q_points +
                 q_point;
        }



        /**
         * Return the quadrature point.
         */
        DEAL_II_HOST_DEVICE
        typename MatrixFree<dim, Number>::point_type &
        get_quadrature_point(const unsigned int cell,
                             const unsigned int q_point) const
        {
          Assert(precomputed_data[0].n_cells ==
                   precomputed_data[0].q_points.extent(1),
                 ExcInternalError());
          AssertIndexRange(cell, precomputed_data[0].n_cells);
          AssertIndexRange(q_point, n_q_points);
          Assert(n_q_points == precomputed_data[0].q_points.extent(0),
                 ExcInternalError());
          return precomputed_data[0].q_points(q_point, cell);
        }

        /**
         * Return the current cell index corresponding to the given @p batch_cell_index.
         * The index returned is only unique for a given MPI process.
         */
        DEAL_II_HOST_DEVICE unsigned int
        current_cell_index(const int batch_cell_index) const
        {
          AssertIndexRange(batch_cell_index, n_cells_in_current_batch);
          return batch_index * n_cells_per_batch + batch_cell_index;
        }

        /**
         * Apply the given functor to each quadrature point of each cell in
         * this batch, exposing the cell index @p func can use directly with
         * local_q_point_id()/get_quadrature_point() or to index per-cell
         * data (e.g. a coefficient, as in step-37's local_apply()) -- the
         * batch-local translation done by current_cell_index() happens
         * once per iteration internally, so @p func doesn't need to call it
         * itself. This keeps user code close to the CPU
         * MatrixFree::cell_loop() pattern, where a cell_range of actual
         * cell indices is what gets iterated, not a batch-local one.
         *
         * @p func needs to define
         * \code
         * DEAL_II_HOST_DEVICE
         * void operator()(const unsigned int &cell,
         *                 const unsigned int &q_point) const;
         * \endcode
         */
        template <typename Functor>
        DEAL_II_HOST_DEVICE void
        for_each_q_point_in_cell_range(const Functor &func) const
        {
          Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team_member,
                                    n_cells_in_current_batch * n_q_points),
            [&](const int batch_q_index) {
              const int batch_cell_index = batch_q_index / n_q_points;
              const int local_q_point    = batch_q_index % n_q_points;

              func(current_cell_index(batch_cell_index), local_q_point);
            });
          team_member.team_barrier();
        }
      };

      /**
       * This struct performs the evaluation of function values and gradients
       * for tensor-product finite elements. There are two specialized
       * implementation classes FEEvaluationImplCollocation (for Gauss-Lobatto
       * elements where the nodal points and the quadrature points coincide and
       * the 'values' operation is identity) and
       * FEEvaluationImplTransformToCollocation (which can be transformed to a
       * collocation space and can then use the identity in these spaces), which
       * both allow for shorter code.
       */
      template <int dim, int fe_degree, int n_q_points_1d, typename Number>
      struct FEEvaluationImpl
      {
        using SharedView =
          Kokkos::View<Number *,
                       MemorySpace::Default::kokkos_space::execution_space::
                         scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;


        DEAL_II_HOST_DEVICE static void
        evaluate(const unsigned int                     dof_handler_index,
                 const unsigned int                     n_components,
                 const EvaluationFlags::EvaluationFlags evaluation_flag,
                 const BatchedData<dim, Number>        *data)
        {
          if (evaluation_flag == EvaluationFlags::nothing)
            return;

          // the evaluator does not need temporary storage since no in-place
          // operation takes place in this function
          auto scratch_for_eval =
            Kokkos::subview(data->shared_data[dof_handler_index].scratch_pad,
                            Kokkos::make_pair(0, 0));
          EvaluatorTensorProduct<EvaluatorVariant::evaluate_general,
                                 dim,
                                 fe_degree + 1,
                                 n_q_points_1d,
                                 Number>
            eval(
              data->team_member,
              data->shared_data[dof_handler_index].shape_values,
              data->shared_data[dof_handler_index].shape_gradients,
              SharedView(), // co_shape_gradients are not needed for evaluation
              scratch_for_eval,
              data->n_cells_in_current_batch);

          for (unsigned int c = 0; c < n_components; ++c)
            {
              auto u =
                Kokkos::subview(data->shared_data[dof_handler_index].values,
                                Kokkos::ALL,
                                c);
              auto grad_u =
                Kokkos::subview(data->shared_data[dof_handler_index].gradients,
                                Kokkos::ALL,
                                Kokkos::ALL,
                                c);

              if constexpr (dim == 1)
                {
                  auto temp = Kokkos::subview(
                    data->shared_data[dof_handler_index].scratch_pad,
                    Kokkos::make_pair(0, n_q_points_1d));

                  if (evaluation_flag & EvaluationFlags::gradients)
                    eval.template gradients<0, true, false, false>(
                      u, Kokkos::subview(grad_u, Kokkos::ALL, 0));
                  if (evaluation_flag & EvaluationFlags::values)
                    {
                      eval.template values<0, true, false, false>(u, temp);
                      populate_view<false>(data->team_member,
                                           u,
                                           temp,
                                           data->n_cells_in_current_batch *
                                             n_q_points_1d);
                    }
                }
              else if constexpr (dim == 2)
                {
                  const int temp_size = data->n_cells_in_current_batch *
                                        (fe_degree + 1) * n_q_points_1d;
                  auto temp = Kokkos::subview(
                    data->shared_data[dof_handler_index].scratch_pad,
                    Kokkos::make_pair(0, temp_size));

                  // grad x
                  if (evaluation_flag & EvaluationFlags::gradients)
                    {
                      eval.template gradients<0, true, false, false>(u, temp);
                      eval.template values<1, true, false, false>(
                        temp, Kokkos::subview(grad_u, Kokkos::ALL, 0));
                    }

                  // grad y
                  eval.template values<0, true, false, false>(u, temp);
                  if (evaluation_flag & EvaluationFlags::gradients)
                    eval.template gradients<1, true, false, false>(
                      temp, Kokkos::subview(grad_u, Kokkos::ALL, 1));

                  // val: can use values applied in x
                  if (evaluation_flag & EvaluationFlags::values)
                    eval.template values<1, true, false, false>(temp, u);
                }
              else if constexpr (dim == 3)
                {
                  const int temp1_size = data->n_cells_in_current_batch *
                                         Utilities::pow(fe_degree + 1, 2) *
                                         n_q_points_1d;
                  const int temp2_size = data->n_cells_in_current_batch *
                                         Utilities::pow(n_q_points_1d, 2) *
                                         (fe_degree + 1);

                  auto temp1 = Kokkos::subview(
                    data->shared_data[dof_handler_index].scratch_pad,
                    Kokkos::make_pair(0, temp1_size));
                  auto temp2 = Kokkos::subview(
                    data->shared_data[dof_handler_index].scratch_pad,
                    Kokkos::make_pair(temp1_size, temp1_size + temp2_size));

                  if (evaluation_flag & EvaluationFlags::gradients)
                    {
                      // grad x
                      eval.template gradients<0, true, false, false>(u, temp1);
                      eval.template values<1, true, false, false>(temp1, temp2);
                      eval.template values<2, true, false, false>(
                        temp2, Kokkos::subview(grad_u, Kokkos::ALL, 0));
                    }

                  // grad y
                  eval.template values<0, true, false, false>(u, temp1);
                  if (evaluation_flag & EvaluationFlags::gradients)
                    {
                      eval.template gradients<1, true, false, false>(temp1,
                                                                     temp2);
                      eval.template values<2, true, false, false>(
                        temp2, Kokkos::subview(grad_u, Kokkos::ALL, 1));
                    }

                  // grad z: can use the values applied in x direction stored
                  // in temp1
                  eval.template values<1, true, false, false>(temp1, temp2);
                  if (evaluation_flag & EvaluationFlags::gradients)
                    eval.template gradients<2, true, false, false>(
                      temp2, Kokkos::subview(grad_u, Kokkos::ALL, 2));

                  // val: can use the values applied in x & y direction stored
                  // in temp2
                  if (evaluation_flag & EvaluationFlags::values)
                    eval.template values<2, true, false, false>(temp2, u);
                }
              else
                Assert(false, ExcMessage("dim must not exceed 3!"));
            }
        }



        DEAL_II_HOST_DEVICE static void
        integrate(const unsigned int                     dof_handler_index,
                  const unsigned int                     n_components,
                  const EvaluationFlags::EvaluationFlags integration_flag,
                  const BatchedData<dim, Number>        *data)
        {
          if (integration_flag == EvaluationFlags::nothing)
            return;

          // the evaluator does not need temporary storage since no in-place
          // operation takes place in this function
          auto scratch_for_eval =
            Kokkos::subview(data->shared_data[dof_handler_index].scratch_pad,
                            Kokkos::make_pair(0, 0));
          EvaluatorTensorProduct<EvaluatorVariant::evaluate_general,
                                 dim,
                                 fe_degree + 1,
                                 n_q_points_1d,
                                 Number>
            eval(
              data->team_member,
              data->shared_data[dof_handler_index].shape_values,
              data->shared_data[dof_handler_index].shape_gradients,
              SharedView(), // co_shape_gradients are not needed for evaluation
              scratch_for_eval,
              data->n_cells_in_current_batch);

          for (unsigned int c = 0; c < n_components; ++c)
            {
              auto u =
                Kokkos::subview(data->shared_data[dof_handler_index].values,
                                Kokkos::ALL,
                                c);
              auto grad_u =
                Kokkos::subview(data->shared_data[dof_handler_index].gradients,
                                Kokkos::ALL,
                                Kokkos::ALL,
                                c);

              if constexpr (dim == 1)
                {
                  auto temp = Kokkos::subview(
                    data->shared_data[dof_handler_index].scratch_pad,
                    Kokkos::make_pair(0,
                                      data->n_cells_in_current_batch *
                                        (fe_degree + 1)));

                  if ((integration_flag & EvaluationFlags::values) &&
                      !(integration_flag & EvaluationFlags::gradients))
                    {
                      eval.template values<0, false, false, false>(u, temp);
                      populate_view<false>(data->team_member,
                                           u,
                                           temp,
                                           data->n_cells_in_current_batch *
                                             (fe_degree + 1));
                    }
                  if (integration_flag & EvaluationFlags::gradients)
                    {
                      if (integration_flag & EvaluationFlags::values)
                        {
                          eval.template values<0, false, false, false>(u, temp);
                          eval.template gradients<0, false, true, false>(
                            Kokkos::subview(grad_u, Kokkos::ALL, 0), temp);
                          populate_view<false>(data->team_member,
                                               u,
                                               temp,
                                               data->n_cells_in_current_batch *
                                                 (fe_degree + 1));
                        }
                      else
                        eval.template gradients<0, false, false, false>(
                          Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
                    }
                }
              else if constexpr (dim == 2)
                {
                  const int temp_size = data->n_cells_in_current_batch *
                                        (fe_degree + 1) * n_q_points_1d;
                  auto temp = Kokkos::subview(
                    data->shared_data[dof_handler_index].scratch_pad,
                    Kokkos::make_pair(0, temp_size));

                  if ((integration_flag & EvaluationFlags::values) &&
                      !(integration_flag & EvaluationFlags::gradients))
                    {
                      eval.template values<1, false, false, false>(u, temp);
                      eval.template values<0, false, false, false>(temp, u);
                    }
                  if (integration_flag & EvaluationFlags::gradients)
                    {
                      eval.template gradients<1, false, false, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 1), temp);
                      if (integration_flag & EvaluationFlags::values)
                        eval.template values<1, false, true, false>(u, temp);
                      eval.template values<0, false, false, false>(temp, u);
                      eval.template values<1, false, false, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 0), temp);
                      eval.template gradients<0, false, true, false>(temp, u);
                    }
                }
              else if constexpr (dim == 3)
                {
                  const int temp1_size = data->n_cells_in_current_batch *
                                         Utilities::pow(n_q_points_1d, 2) *
                                         (fe_degree + 1);
                  const int temp2_size = data->n_cells_in_current_batch *
                                         Utilities::pow(fe_degree + 1, 2) *
                                         n_q_points_1d;

                  auto temp1 = Kokkos::subview(
                    data->shared_data[dof_handler_index].scratch_pad,
                    Kokkos::make_pair(0, temp1_size));
                  auto temp2 = Kokkos::subview(
                    data->shared_data[dof_handler_index].scratch_pad,
                    Kokkos::make_pair(temp1_size, temp1_size + temp2_size));

                  if ((integration_flag & EvaluationFlags::values) &&
                      !(integration_flag & EvaluationFlags::gradients))
                    {
                      eval.template values<2, false, false, false>(u, temp1);
                      eval.template values<1, false, false, false>(temp1,
                                                                   temp2);
                      eval.template values<0, false, false, false>(temp2, u);
                    }
                  if (integration_flag & EvaluationFlags::gradients)
                    {
                      eval.template gradients<2, false, false, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 2), temp1);
                      if (integration_flag & EvaluationFlags::values)
                        eval.template values<2, false, true, false>(u, temp1);
                      eval.template values<1, false, false, false>(temp1,
                                                                   temp2);
                      eval.template values<2, false, false, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 1), temp1);
                      eval.template gradients<1, false, true, false>(temp1,
                                                                     temp2);
                      eval.template values<0, false, false, false>(temp2, u);
                      eval.template values<2, false, false, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 0), temp1);
                      eval.template values<1, false, false, false>(temp1,
                                                                   temp2);
                      eval.template gradients<0, false, true, false>(temp2, u);
                    }
                }
              else
                Assert(false, ExcMessage("dim must not exceed 3!"));
            }
        }
      };


      /**
       * This struct performs the evaluation of function values and gradients
       * for tensor-product finite elements. This is a specialization for
       * elements where the nodal points coincide with the quadrature points
       * like FE_Q shape functions on Gauss-Lobatto elements integrated with
       * Gauss-Lobatto quadrature. The assumption of this class is that the
       * shape 'values' operation is identity, which allows us to write shorter
       * code.
       *
       * In literature, this form of evaluation is often called spectral
       * evaluation, spectral collocation or simply collocation, meaning the
       * same location for shape functions and evaluation space (quadrature
       * points).
       */
      template <int dim, int fe_degree, typename Number>
      struct FEEvaluationImplCollocation
      {
        using SharedView =
          Kokkos::View<Number *,
                       MemorySpace::Default::kokkos_space::execution_space::
                         scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

        DEAL_II_HOST_DEVICE static void
        evaluate(const unsigned int                     dof_handler_index,
                 const unsigned int                     n_components,
                 const EvaluationFlags::EvaluationFlags evaluation_flag,
                 const BatchedData<dim, Number>        *data)
        {
          Assert(dim >= 1 && dim <= 3, ExcMessage("dim must be 1, 2, or 3"));

          // since the dof values have already been stored in
          // shared_data->values, there is nothing to do if the gradients are
          // not required
          if (!(evaluation_flag & EvaluationFlags::gradients))
            return;

          // the evaluator does not need temporary storage since no in-place
          // operation takes place in this function
          auto scratch_for_eval =
            Kokkos::subview(data->shared_data[dof_handler_index].scratch_pad,
                            Kokkos::make_pair(0, 0));

          EvaluatorTensorProduct<EvaluatorVariant::evaluate_general,
                                 dim,
                                 fe_degree + 1,
                                 fe_degree + 1,
                                 Number>
            eval(data->team_member,
                 SharedView(), // shape_values are identity
                 SharedView(), // shape_gradients are not used
                 data->shared_data[dof_handler_index].co_shape_gradients,
                 scratch_for_eval,
                 data->n_cells_in_current_batch);

          for (unsigned int c = 0; c < n_components; ++c)
            {
              auto u =
                Kokkos::subview(data->shared_data[dof_handler_index].values,
                                Kokkos::ALL,
                                c);
              auto grad_u =
                Kokkos::subview(data->shared_data[dof_handler_index].gradients,
                                Kokkos::ALL,
                                Kokkos::ALL,
                                c);
              // apply derivatives in collocation space fused along all dim
              // directions
              eval.template co_gradients<false, false>(u, grad_u);
            }
        }


        DEAL_II_HOST_DEVICE static void
        integrate(const unsigned int                     dof_handler_index,
                  const unsigned int                     n_components,
                  const EvaluationFlags::EvaluationFlags integration_flag,
                  const BatchedData<dim, Number>        *data)
        {
          Assert(dim >= 1 && dim <= 3, ExcMessage("dim must be 1, 2, or 3"));

          // since the quad values have already been stored in
          // shared_data->values, there is nothing to do if the gradients are
          // not required
          if (!(integration_flag & EvaluationFlags::gradients))
            return;

          // the evaluator does not need temporary storage since no in-place
          // operation takes place in this function
          auto scratch_for_eval =
            Kokkos::subview(data->shared_data[dof_handler_index].scratch_pad,
                            Kokkos::make_pair(0, 0));

          EvaluatorTensorProduct<EvaluatorVariant::evaluate_general,
                                 dim,
                                 fe_degree + 1,
                                 fe_degree + 1,
                                 Number>
            eval(data->team_member,
                 SharedView(), // shape_values are identity
                 SharedView(), // shape_gradients are not used
                 data->shared_data[dof_handler_index].co_shape_gradients,
                 scratch_for_eval,
                 data->n_cells_in_current_batch);

          for (unsigned int c = 0; c < n_components; ++c)
            {
              auto u =
                Kokkos::subview(data->shared_data[dof_handler_index].values,
                                Kokkos::ALL,
                                c);
              auto grad_u =
                Kokkos::subview(data->shared_data[dof_handler_index].gradients,
                                Kokkos::ALL,
                                Kokkos::ALL,
                                c);

              // apply derivatives in collocation space fused along all dim
              // directions
              if (integration_flag & EvaluationFlags::values)
                eval.template co_gradients<true, true>(grad_u, u);
              else
                eval.template co_gradients<true, false>(grad_u, u);
            }
        }
      };

      /**
       * This struct performs the evaluation of function values and gradients
       * for tensor-product finite elements. This is a specialization for
       * symmetric basis functions about the mid point 0.5 of the unit
       * interval with the same number of quadrature points as degrees of
       * freedom. In that case, we can first transform the basis to one that has
       * the nodal points in the quadrature points (i.e., the collocation space)
       * and then perform the evaluation of the first and second derivatives in
       * this transformed space, using the identity operation for the shape
       * values.
       */
      template <int dim, int fe_degree, int n_q_points_1d, typename Number>
      struct FEEvaluationImplTransformToCollocation
      {
        using SharedView =
          Kokkos::View<Number *,
                       MemorySpace::Default::kokkos_space::execution_space::
                         scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

        DEAL_II_HOST_DEVICE static void
        evaluate(const unsigned int                     dof_handler_index,
                 const unsigned int                     n_components,
                 const EvaluationFlags::EvaluationFlags evaluation_flag,
                 const BatchedData<dim, Number>        *data)
        {
          // the evaluator does not need temporary storage since no in-place
          // operation takes place in this function as each in-place operation
          // requires populate_view() copy and extra team_barrier(). We call
          // populate_view() manually only when necessary, and this bring around
          // 5% performance improvement in 3D
          auto scratch_for_eval =
            Kokkos::subview(data->shared_data[dof_handler_index].scratch_pad,
                            Kokkos::make_pair(0, 0));

          EvaluatorTensorProduct<EvaluatorVariant::evaluate_general,
                                 dim,
                                 fe_degree + 1,
                                 n_q_points_1d,
                                 Number>
            eval(data->team_member,
                 data->shared_data[dof_handler_index].shape_values,
                 SharedView(), // shape_gradients are not used
                 data->shared_data[dof_handler_index].co_shape_gradients,
                 scratch_for_eval,
                 data->n_cells_in_current_batch);

          for (unsigned int c = 0; c < n_components; ++c)
            {
              auto u =
                Kokkos::subview(data->shared_data[dof_handler_index].values,
                                Kokkos::ALL,
                                c);
              auto grad_u =
                Kokkos::subview(data->shared_data[dof_handler_index].gradients,
                                Kokkos::ALL,
                                Kokkos::ALL,
                                c);

              if constexpr (dim == 1)
                {
                  const auto &temp =
                    data->shared_data[dof_handler_index].scratch_pad;

                  eval.template values<0, true, false, false>(u, temp);

                  populate_view<false>(data->team_member,
                                       u,
                                       temp,
                                       data->n_cells_in_current_batch *
                                         data->n_q_points);
                }
              else if constexpr (dim == 2)
                {
                  const auto &temp =
                    data->shared_data[dof_handler_index].scratch_pad;

                  eval.template values<0, true, false, false>(u, temp);
                  eval.template values<1, true, false, false>(temp, u);
                }
              else // dim == 3
                {
                  const auto &temp =
                    data->shared_data[dof_handler_index].scratch_pad;

                  eval.template values<0, true, false, false>(u, temp);
                  eval.template values<1, true, false, false>(temp, u);
                  eval.template values<2, true, false, false>(u, temp);

                  populate_view<false>(data->team_member,
                                       u,
                                       temp,
                                       data->n_cells_in_current_batch *
                                         data->n_q_points);
                }

              // apply derivatives in collocation space fused along all dim
              // directions
              if (evaluation_flag & EvaluationFlags::gradients)
                eval.template co_gradients<false, false>(u, grad_u);
            }
        }


        DEAL_II_HOST_DEVICE static void
        integrate(const unsigned int                     dof_handler_index,
                  const unsigned int                     n_components,
                  const EvaluationFlags::EvaluationFlags integration_flag,
                  const BatchedData<dim, Number>        *data)
        {
          // the evaluator does not need temporary storage since no in-place
          // operation takes place in this function as each in-place operation
          // requires populate_view() copy and extra team_barrier(). We call
          // populate_view() manually only when necessary, and this bring around
          // 5% performance improvement in 3D
          auto scratch_for_eval =
            Kokkos::subview(data->shared_data[dof_handler_index].scratch_pad,
                            Kokkos::make_pair(0, 0));

          EvaluatorTensorProduct<EvaluatorVariant::evaluate_general,
                                 dim,
                                 fe_degree + 1,
                                 n_q_points_1d,
                                 Number>
            eval(data->team_member,
                 data->shared_data[dof_handler_index].shape_values,
                 SharedView(), // shape_gradients are not used
                 data->shared_data[dof_handler_index].co_shape_gradients,
                 scratch_for_eval,
                 data->n_cells_in_current_batch);

          for (unsigned int c = 0; c < n_components; ++c)
            {
              auto u =
                Kokkos::subview(data->shared_data[dof_handler_index].values,
                                Kokkos::ALL,
                                c);
              auto grad_u =
                Kokkos::subview(data->shared_data[dof_handler_index].gradients,
                                Kokkos::ALL,
                                Kokkos::ALL,
                                c);

              // apply derivatives in collocation space fused along all dim
              // directions
              if (integration_flag & EvaluationFlags::gradients)
                {
                  if (integration_flag & EvaluationFlags::values)
                    eval.template co_gradients<true, true>(grad_u, u);
                  else
                    eval.template co_gradients<true, false>(grad_u, u);
                }

              if constexpr (dim == 1)
                {
                  const auto &temp =
                    data->shared_data[dof_handler_index].scratch_pad;

                  eval.template values<0, false, false, false>(u, temp);

                  populate_view<false>(data->team_member,
                                       u,
                                       temp,
                                       data->n_cells_in_current_batch *
                                         (fe_degree + 1));
                }
              else if constexpr (dim == 2)
                {
                  const auto &temp =
                    data->shared_data[dof_handler_index].scratch_pad;

                  eval.template values<1, false, false, false>(u, temp);
                  eval.template values<0, false, false, false>(temp, u);
                }
              else if constexpr (dim == 3)
                {
                  const auto &temp =
                    data->shared_data[dof_handler_index].scratch_pad;

                  eval.template values<2, false, false, false>(u, temp);
                  eval.template values<1, false, false, false>(temp, u);
                  eval.template values<0, false, false, false>(u, temp);

                  populate_view<false>(data->team_member,
                                       u,
                                       temp,
                                       data->n_cells_in_current_batch *
                                         Utilities::pow(fe_degree + 1, dim));
                }
              else
                Assert(false, ExcMessage("dim must not exceed 3!"));
            }
        }
      };

    } // namespace batched
  }   // end of namespace internal
} // end of namespace Portable


DEAL_II_NAMESPACE_CLOSE

#endif
