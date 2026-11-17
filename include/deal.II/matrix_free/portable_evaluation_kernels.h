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
         * Allocate Kokkos Views for a batch of @p batch_size cells and copy
         * the shape data from the precomputed data to shared memory. The
         * values, gradients and scratch pad of the cells are stored cell
         * after cell.
         */
        DEAL_II_HOST_DEVICE
        void
        reinit(
          const TeamHandle  &team_handle,
          const unsigned int n_q_points_per_cell,
          const int          batch_size,
          const typename MatrixFree<dim, Number>::PrecomputedData &gpu_data)
        {
          const unsigned int n_points = batch_size * n_q_points_per_cell;

          values    = SharedViewValues(team_handle.team_shmem(),
                                       n_points,
                                       gpu_data.n_components);
          gradients = SharedViewGradients(team_handle.team_shmem(),
                                          n_points,
                                          dim,
                                          gpu_data.n_components);
          scratch_pad =
            SharedViewScratchPad(team_handle.team_shmem(),
                                 batch_size * gpu_data.scratch_pad_size);

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
         * Number of bytes of team scratch memory that reinit() allocates for
         * the same arguments.
         */
        static std::size_t
        shmem_size(
          const unsigned int n_q_points_per_cell,
          const int          batch_size,
          const typename MatrixFree<dim, Number>::PrecomputedData &gpu_data)
        {
          using ElementType =
            ::dealii::internal::MatrixFreeFunctions::ElementType;

          const unsigned int n_points = batch_size * n_q_points_per_cell;

          std::size_t result =
            SharedViewValues::shmem_size(n_points, gpu_data.n_components) +
            SharedViewGradients::shmem_size(n_points,
                                            dim,
                                            gpu_data.n_components) +
            SharedViewScratchPad::shmem_size(batch_size *
                                             gpu_data.scratch_pad_size);

          if (gpu_data.element_type >= ElementType::tensor_symmetric)
            result +=
              SharedViewShapeData::shmem_size(gpu_data.shape_values.size());

          if (gpu_data.element_type <= ElementType::tensor_symmetric)
            result += SharedViewShapeData::shmem_size(
              gpu_data.co_shape_gradients.size());
          else
            result +=
              SharedViewShapeData::shmem_size(gpu_data.shape_gradients.size());

          return result;
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
       * Same as Portable::MatrixFree::Data, but with the above SharedData
       * struct and with a range of cells handled by the team instead of a
       * single cell index.
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

        /**
         * Number of quadrature points per cell.
         */
        const unsigned int n_q_points;

        /**
         * Number of DoFHandler objects, i.e., of used entries in
         * @p precomputed_data and @p shared_data.
         */
        const unsigned int n_dof_handler;

        /**
         * Half-open range [first, second) of the cells handled by this team,
         * in the same index space as Portable::MatrixFree::Data::cell_index.
         * The last batch of a cell loop may be shorter than the others.
         */
        const Kokkos::pair<int, int> cell_range;

        /**
         * Precomputed data of each DoFHandler (shape data, mapping data, dof
         * indices) in global memory.
         */
        const Kokkos::Array<typename MatrixFree<dim, Number>::PrecomputedData,
                            n_max_dof_handlers> &precomputed_data;

        /**
         * Shared memory of each DoFHandler: staged shape data, the dof or
         * quadrature point values and gradients of all cells in @p cell_range
         * stored cell after cell, and scratch space.
         */
        Kokkos::Array<SharedData<dim, Number>, n_max_dof_handlers> &shared_data;

        /**
         * Number of cells in @p cell_range.
         */
        DEAL_II_HOST_DEVICE int
        batch_size() const
        {
          Assert(cell_range.first <= cell_range.second, ExcInternalError());
          return cell_range.second - cell_range.first;
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

          const int batch_size = data->batch_size();

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
              batch_size);

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
                    Kokkos::make_pair(0, batch_size * n_q_points_1d));

                  if (evaluation_flag & EvaluationFlags::gradients)
                    eval.template derivative<0, true, false, false>(
                      u, Kokkos::subview(grad_u, Kokkos::ALL, 0));
                  if (evaluation_flag & EvaluationFlags::values)
                    {
                      eval.template value<0, true, false, false>(u, temp);
                      populate_view<false>(data->team_member,
                                           u,
                                           temp,
                                           batch_size * n_q_points_1d);
                    }
                }
              else if constexpr (dim == 2)
                {
                  const int temp_size =
                    batch_size * (fe_degree + 1) * n_q_points_1d;
                  auto temp = Kokkos::subview(
                    data->shared_data[dof_handler_index].scratch_pad,
                    Kokkos::make_pair(0, temp_size));

                  // grad x
                  if (evaluation_flag & EvaluationFlags::gradients)
                    {
                      eval.template derivative<0, true, false, false>(u, temp);
                      eval.template value<1, true, false, false>(
                        temp, Kokkos::subview(grad_u, Kokkos::ALL, 0));
                    }

                  // grad y
                  eval.template value<0, true, false, false>(u, temp);
                  if (evaluation_flag & EvaluationFlags::gradients)
                    eval.template derivative<1, true, false, false>(
                      temp, Kokkos::subview(grad_u, Kokkos::ALL, 1));

                  // val: can use values applied in x
                  if (evaluation_flag & EvaluationFlags::values)
                    eval.template value<1, true, false, false>(temp, u);
                }
              else if constexpr (dim == 3)
                {
                  const int temp1_size = batch_size *
                                         Utilities::pow(fe_degree + 1, 2) *
                                         n_q_points_1d;
                  const int temp2_size = batch_size *
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
                      eval.template derivative<0, true, false, false>(u, temp1);
                      eval.template value<1, true, false, false>(temp1, temp2);
                      eval.template value<2, true, false, false>(
                        temp2, Kokkos::subview(grad_u, Kokkos::ALL, 0));
                    }

                  // grad y
                  eval.template value<0, true, false, false>(u, temp1);
                  if (evaluation_flag & EvaluationFlags::gradients)
                    {
                      eval.template derivative<1, true, false, false>(temp1,
                                                                      temp2);
                      eval.template value<2, true, false, false>(
                        temp2, Kokkos::subview(grad_u, Kokkos::ALL, 1));
                    }

                  // grad z: can use the values applied in x direction stored
                  // in temp1
                  eval.template value<1, true, false, false>(temp1, temp2);
                  if (evaluation_flag & EvaluationFlags::gradients)
                    eval.template derivative<2, true, false, false>(
                      temp2, Kokkos::subview(grad_u, Kokkos::ALL, 2));

                  // val: can use the values applied in x & y direction stored
                  // in temp2
                  if (evaluation_flag & EvaluationFlags::values)
                    eval.template value<2, true, false, false>(temp2, u);
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

          const int batch_size = data->batch_size();

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
              batch_size);

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
                    Kokkos::make_pair(0, batch_size * (fe_degree + 1)));

                  if ((integration_flag & EvaluationFlags::values) &&
                      !(integration_flag & EvaluationFlags::gradients))
                    {
                      eval.template value<0, false, false, false>(u, temp);
                      populate_view<false>(data->team_member,
                                           u,
                                           temp,
                                           batch_size * (fe_degree + 1));
                    }
                  if (integration_flag & EvaluationFlags::gradients)
                    {
                      if (integration_flag & EvaluationFlags::values)
                        {
                          eval.template value<0, false, false, false>(u, temp);
                          eval.template derivative<0, false, true, false>(
                            Kokkos::subview(grad_u, Kokkos::ALL, 0), temp);
                          populate_view<false>(data->team_member,
                                               u,
                                               temp,
                                               batch_size * (fe_degree + 1));
                        }
                      else
                        eval.template derivative<0, false, false, false>(
                          Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
                    }
                }
              else if constexpr (dim == 2)
                {
                  const int temp_size =
                    batch_size * (fe_degree + 1) * n_q_points_1d;
                  auto temp = Kokkos::subview(
                    data->shared_data[dof_handler_index].scratch_pad,
                    Kokkos::make_pair(0, temp_size));

                  if ((integration_flag & EvaluationFlags::values) &&
                      !(integration_flag & EvaluationFlags::gradients))
                    {
                      eval.template value<1, false, false, false>(u, temp);
                      eval.template value<0, false, false, false>(temp, u);
                    }
                  if (integration_flag & EvaluationFlags::gradients)
                    {
                      eval.template derivative<1, false, false, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 1), temp);
                      if (integration_flag & EvaluationFlags::values)
                        eval.template value<1, false, true, false>(u, temp);
                      eval.template value<0, false, false, false>(temp, u);
                      eval.template value<1, false, false, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 0), temp);
                      eval.template derivative<0, false, true, false>(temp, u);
                    }
                }
              else if constexpr (dim == 3)
                {
                  const int temp1_size = batch_size *
                                         Utilities::pow(n_q_points_1d, 2) *
                                         (fe_degree + 1);
                  const int temp2_size = batch_size *
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
                      eval.template value<2, false, false, false>(u, temp1);
                      eval.template value<1, false, false, false>(temp1, temp2);
                      eval.template value<0, false, false, false>(temp2, u);
                    }
                  if (integration_flag & EvaluationFlags::gradients)
                    {
                      eval.template derivative<2, false, false, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 2), temp1);
                      if (integration_flag & EvaluationFlags::values)
                        eval.template value<2, false, true, false>(u, temp1);
                      eval.template value<1, false, false, false>(temp1, temp2);
                      eval.template value<2, false, false, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 1), temp1);
                      eval.template derivative<1, false, true, false>(temp1,
                                                                      temp2);
                      eval.template value<0, false, false, false>(temp2, u);
                      eval.template value<2, false, false, false>(
                        Kokkos::subview(grad_u, Kokkos::ALL, 0), temp1);
                      eval.template value<1, false, false, false>(temp1, temp2);
                      eval.template derivative<0, false, true, false>(temp2, u);
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

          const int batch_size = data->batch_size();

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
                 batch_size);

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
              eval.template co_gradient<false, false>(u, grad_u);
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

          const int batch_size = data->batch_size();

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
                 batch_size);

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
                eval.template co_gradient<true, true>(grad_u, u);
              else
                eval.template co_gradient<true, false>(grad_u, u);
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
          const int batch_size = data->batch_size();

          // the evaluator does not need temporary storage since no in-place
          // operation takes place in this function as each in-place operation
          // requires populate_view() copy and extra team_barrier(). We call
          // populate_view() manually only when necessary, and this bring around
          // 5% performance improvement
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
                 batch_size);

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

                  eval.template value<0, true, false, false>(u, temp);

                  populate_view<false>(data->team_member,
                                       u,
                                       temp,
                                       batch_size * data->n_q_points);
                }
              else if constexpr (dim == 2)
                {
                  const auto &temp =
                    data->shared_data[dof_handler_index].scratch_pad;

                  eval.template value<0, true, false, false>(u, temp);
                  eval.template value<1, true, false, false>(temp, u);
                }
              else // dim == 3
                {
                  const auto &temp =
                    data->shared_data[dof_handler_index].scratch_pad;

                  eval.template value<0, true, false, false>(u, temp);
                  eval.template value<1, true, false, false>(temp, u);
                  eval.template value<2, true, false, false>(u, temp);

                  populate_view<false>(data->team_member,
                                       u,
                                       temp,
                                       batch_size * data->n_q_points);
                }

              // apply derivatives in collocation space fused along all dim
              // directions
              if (evaluation_flag & EvaluationFlags::gradients)
                eval.template co_gradient<false, false>(u, grad_u);
            }
        }


        DEAL_II_HOST_DEVICE static void
        integrate(const unsigned int                     dof_handler_index,
                  const unsigned int                     n_components,
                  const EvaluationFlags::EvaluationFlags integration_flag,
                  const BatchedData<dim, Number>        *data)
        {
          const int batch_size = data->batch_size();

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
                 batch_size);

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
                    eval.template co_gradient<true, true>(grad_u, u);
                  else
                    eval.template co_gradient<true, false>(grad_u, u);
                }

              if constexpr (dim == 1)
                {
                  const auto &temp =
                    data->shared_data[dof_handler_index].scratch_pad;

                  eval.template value<0, false, false, false>(u, temp);

                  populate_view<false>(data->team_member,
                                       u,
                                       temp,
                                       batch_size * (fe_degree + 1));
                }
              else if constexpr (dim == 2)
                {
                  const auto &temp =
                    data->shared_data[dof_handler_index].scratch_pad;

                  eval.template value<1, false, false, false>(u, temp);
                  eval.template value<0, false, false, false>(temp, u);
                }
              else if constexpr (dim == 3)
                {
                  const auto &temp =
                    data->shared_data[dof_handler_index].scratch_pad;

                  eval.template value<2, false, false, false>(u, temp);
                  eval.template value<1, false, false, false>(temp, u);
                  eval.template value<0, false, false, false>(u, temp);

                  populate_view<false>(data->team_member,
                                       u,
                                       temp,
                                       batch_size *
                                         Utilities::pow(fe_degree + 1, dim));
                }
              else
                Assert(false, ExcMessage("dim must not exceed 3!"));
            }
        }
      };

    } // namespace batched
  } // end of namespace internal
} // end of namespace Portable


DEAL_II_NAMESPACE_CLOSE

#endif
