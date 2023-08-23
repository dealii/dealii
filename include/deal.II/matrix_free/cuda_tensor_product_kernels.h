// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2023 by the deal.II authors
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


#ifndef dealii_cuda_tensor_product_kernels_h
#define dealii_cuda_tensor_product_kernels_h

#include <deal.II/base/config.h>

#include <deal.II/base/utilities.h>

#include <deal.II/matrix_free/cuda_matrix_free.templates.h>

DEAL_II_NAMESPACE_OPEN


namespace CUDAWrappers
{
  namespace internal
  {
    /**
     * In this namespace, the evaluator routines that evaluate the tensor
     * products are implemented.
     *
     * @ingroup CUDAWrappers
     */
    // TODO: for now only the general variant is implemented
    enum EvaluatorVariant
    {
      evaluate_general,
      evaluate_symmetric,
      evaluate_evenodd
    };



#if KOKKOS_VERSION >= 40000
    /**
     * Helper function for values() and gradients() in 1D
     */
    template <int n_q_points_1d,
              typename Number,
              int  direction,
              bool dof_to_quad,
              bool add,
              bool in_place,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    apply_1d(const Kokkos::TeamPolicy<
               MemorySpace::Default::kokkos_space::execution_space>::member_type
               &team_member,
             const Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
                              shape_data,
             const ViewTypeIn in,
             ViewTypeOut      out)
    {
      Number t[n_q_points_1d];
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_q_points_1d),
                           [&](const int &q) {
                             t[q] = 0;
                             // This loop simply multiplies the shape function
                             // at the quadrature point by the value finite
                             // element coefficient.
                             // FIXME check why using parallel_reduce
                             // ThreadVector is slower
                             for (int k = 0; k < n_q_points_1d; ++k)
                               {
                                 const unsigned int shape_idx =
                                   dof_to_quad ? (q + k * n_q_points_1d) :
                                                 (k + q * n_q_points_1d);
                                 const unsigned int source_idx = k;
                                 t[q] += shape_data[shape_idx] * in(source_idx);
                               }
                           });

      if constexpr (in_place)
        team_member.team_barrier();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_q_points_1d),
                           [&](const int &q) {
                             const unsigned int destination_idx = q;
                             if constexpr (add)
                               Kokkos::atomic_add(&out(destination_idx), t[q]);
                             else
                               out(destination_idx) = t[q];
                           });
    }



    /**
     * Helper function for values() and gradients() in 2D
     */
    template <int n_q_points_1d,
              typename Number,
              int  direction,
              bool dof_to_quad,
              bool add,
              bool in_place,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    apply_2d(const Kokkos::TeamPolicy<
               MemorySpace::Default::kokkos_space::execution_space>::member_type
               &team_member,
             const Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
                              shape_data,
             const ViewTypeIn in,
             ViewTypeOut      out)
    {
      using TeamType = Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type;
      constexpr unsigned int n_q_points = Utilities::pow(n_q_points_1d, 2);

      Number t[n_q_points];
      auto   thread_policy =
        Kokkos::TeamThreadMDRange<Kokkos::Rank<2>, TeamType>(team_member,
                                                             n_q_points_1d,
                                                             n_q_points_1d);
      Kokkos::parallel_for(thread_policy, [&](const int i, const int j) {
        int q_point = i + j * n_q_points_1d;
        t[q_point]  = 0;

        // This loop simply multiplies the shape function at the quadrature
        // point by the value finite element coefficient.
        // FIXME check why using parallel_reduce ThreadVector is slower
        for (int k = 0; k < n_q_points_1d; ++k)
          {
            const unsigned int shape_idx =
              dof_to_quad ? (j + k * n_q_points_1d) : (k + j * n_q_points_1d);
            const unsigned int source_idx = (direction == 0) ?
                                              (k + n_q_points_1d * i) :
                                              (i + n_q_points_1d * k);
            t[q_point] += shape_data[shape_idx] * in(source_idx);
          }
      });

      if (in_place)
        team_member.team_barrier();

      Kokkos::parallel_for(thread_policy, [&](const int i, const int j) {
        const int          q_point = i + j * n_q_points_1d;
        const unsigned int destination_idx =
          (direction == 0) ? (j + n_q_points_1d * i) : (i + n_q_points_1d * j);

        if (add)
          Kokkos::atomic_add(&out(destination_idx), t[q_point]);
        else
          out(destination_idx) = t[q_point];
      });
    }



    /**
     * Helper function for values() and gradients() in 3D
     */
    template <int n_q_points_1d,
              typename Number,
              int  direction,
              bool dof_to_quad,
              bool add,
              bool in_place,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    apply_3d(const Kokkos::TeamPolicy<
               MemorySpace::Default::kokkos_space::execution_space>::member_type
               &team_member,
             const Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
                              shape_data,
             const ViewTypeIn in,
             ViewTypeOut      out)
    {
      using TeamType = Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type;
      constexpr unsigned int n_q_points = Utilities::pow(n_q_points_1d, 3);

      Number t[n_q_points];
      auto thread_policy = Kokkos::TeamThreadMDRange<Kokkos::Rank<3>, TeamType>(
        team_member, n_q_points_1d, n_q_points_1d, n_q_points_1d);
      Kokkos::parallel_for(
        thread_policy, [&](const int i, const int j, const int q) {
          const int q_point =
            i + j * n_q_points_1d + q * n_q_points_1d * n_q_points_1d;
          t[q_point] = 0;

          // This loop simply multiplies the shape function at the quadrature
          // point by the value finite element coefficient.
          // FIXME check why using parallel_reduce ThreadVector is slower
          for (int k = 0; k < n_q_points_1d; ++k)
            {
              const unsigned int shape_idx =
                dof_to_quad ? (q + k * n_q_points_1d) : (k + q * n_q_points_1d);
              const unsigned int source_idx =
                (direction == 0) ?
                  (k + n_q_points_1d * (i + n_q_points_1d * j)) :
                (direction == 1) ?
                  (i + n_q_points_1d * (k + n_q_points_1d * j)) :
                  (i + n_q_points_1d * (j + n_q_points_1d * k));
              t[q_point] += shape_data[shape_idx] * in(source_idx);
            }
        });

      if (in_place)
        team_member.team_barrier();

      Kokkos::parallel_for(
        thread_policy, [&](const int i, const int j, const int q) {
          const int q_point =
            i + j * n_q_points_1d + q * n_q_points_1d * n_q_points_1d;
          const unsigned int destination_idx =
            (direction == 0) ? (q + n_q_points_1d * (i + n_q_points_1d * j)) :
            (direction == 1) ? (i + n_q_points_1d * (q + n_q_points_1d * j)) :
                               (i + n_q_points_1d * (j + n_q_points_1d * q));

          if (add)
            Kokkos::atomic_add(&out(destination_idx), t[q_point]);
          else
            out(destination_idx) = t[q_point];
        });
    }
#endif



    /**
     * Helper function for values() and gradients().
     */
    template <int dim,
              int n_q_points_1d,
              typename Number,
              int  direction,
              bool dof_to_quad,
              bool add,
              bool in_place,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    apply(const Kokkos::TeamPolicy<
            MemorySpace::Default::kokkos_space::execution_space>::member_type
            &team_member,
          const Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
                           shape_data,
          const ViewTypeIn in,
          ViewTypeOut      out)
    {
#if KOKKOS_VERSION >= 40000
      if constexpr (dim == 1)
        apply_1d<n_q_points_1d, Number, direction, dof_to_quad, add, in_place>(
          team_member, shape_data, in, out);
      if constexpr (dim == 2)
        apply_2d<n_q_points_1d, Number, direction, dof_to_quad, add, in_place>(
          team_member, shape_data, in, out);
      if constexpr (dim == 3)
        apply_3d<n_q_points_1d, Number, direction, dof_to_quad, add, in_place>(
          team_member, shape_data, in, out);
#else
      constexpr unsigned int n_q_points = Utilities::pow(n_q_points_1d, dim);

      Number t[n_q_points];
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team_member, n_q_points),
        [&](const int &q_point) {
          const unsigned int i = (dim == 1) ? 0 : q_point % n_q_points_1d;
          const unsigned int j =
            (dim == 3) ? (q_point / n_q_points_1d) % n_q_points_1d : 0;
          const unsigned int q =
            (dim == 1) ? q_point :
            (dim == 2) ? (q_point / n_q_points_1d) % n_q_points_1d :
                         q_point / (n_q_points_1d * n_q_points_1d);

          // This loop simply multiplies the shape function at the quadrature
          // point by the value finite element coefficient.
          t[q_point] = 0;
          for (int k = 0; k < n_q_points_1d; ++k)
            {
              const unsigned int shape_idx =
                dof_to_quad ? (q + k * n_q_points_1d) : (k + q * n_q_points_1d);
              const unsigned int source_idx =
                (direction == 0) ?
                  (k + n_q_points_1d * (i + n_q_points_1d * j)) :
                (direction == 1) ?
                  (i + n_q_points_1d * (k + n_q_points_1d * j)) :
                  (i + n_q_points_1d * (j + n_q_points_1d * k));
              t[q_point] += shape_data[shape_idx] *
                            (in_place ? out(source_idx) : in(source_idx));
            }
        });

      if (in_place)
        team_member.team_barrier();

      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team_member, n_q_points),
        [&](const int &q_point) {
          const unsigned int i = (dim == 1) ? 0 : q_point % n_q_points_1d;
          const unsigned int j =
            (dim == 3) ? (q_point / n_q_points_1d) % n_q_points_1d : 0;
          const unsigned int q =
            (dim == 1) ? q_point :
            (dim == 2) ? (q_point / n_q_points_1d) % n_q_points_1d :
                         q_point / (n_q_points_1d * n_q_points_1d);

          const unsigned int destination_idx =
            (direction == 0) ? (q + n_q_points_1d * (i + n_q_points_1d * j)) :
            (direction == 1) ? (i + n_q_points_1d * (q + n_q_points_1d * j)) :
                               (i + n_q_points_1d * (j + n_q_points_1d * q));

          if (add)
            Kokkos::atomic_add(&out(destination_idx), t[q_point]);
          else
            out(destination_idx) = t[q_point];
        });
#endif
    }


    /**
     * Generic evaluator framework.
     *
     * @ingroup CUDAWrappers
     */
    template <EvaluatorVariant variant,
              int              dim,
              int              fe_degree,
              int              n_q_points_1d,
              typename Number>
    struct EvaluatorTensorProduct
    {};



    /**
     * Internal evaluator for 1d-3d shape function using the tensor product form
     * of the basis functions.
     *
     * @ingroup CUDAWrappers
     */
    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    struct EvaluatorTensorProduct<evaluate_general,
                                  dim,
                                  fe_degree,
                                  n_q_points_1d,
                                  Number>
    {
      using TeamHandle = Kokkos::TeamPolicy<
        MemorySpace::Default::kokkos_space::execution_space>::member_type;

      DEAL_II_HOST_DEVICE
      EvaluatorTensorProduct(
        const TeamHandle                                          &team_member,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space> shape_values,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
          shape_gradients,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
          co_shape_gradients);

      /**
       * Evaluate the values of a finite element function at the quadrature
       * points.
       */
      template <int  direction,
                bool dof_to_quad,
                bool add,
                bool in_place,
                typename ViewTypeIn,
                typename ViewTypeOut>
      DEAL_II_HOST_DEVICE void
      values(const ViewTypeIn in, ViewTypeOut out) const;

      /**
       * Evaluate the gradient of a finite element function at the quadrature
       * points for a given @p direction.
       */
      template <int  direction,
                bool dof_to_quad,
                bool add,
                bool in_place,
                typename ViewTypeIn,
                typename ViewTypeOut>
      DEAL_II_HOST_DEVICE void
      gradients(const ViewTypeIn in, ViewTypeOut out) const;

      /**
       * Evaluate the gradient of a finite element function at the quadrature
       * points for a given @p direction for collocation methods.
       */
      template <int  direction,
                bool dof_to_quad,
                bool add,
                bool in_place,
                typename ViewTypeIn,
                typename ViewTypeOut>
      DEAL_II_HOST_DEVICE void
      co_gradients(const ViewTypeIn in, ViewTypeOut out) const;

      /**
       * Evaluate the finite element function at the quadrature points.
       */
      template <typename ViewType>
      DEAL_II_HOST_DEVICE void
      value_at_quad_pts(ViewType u);

      /**
       * Helper function for integrate(). Integrate the finite element function.
       */
      template <typename ViewType>
      DEAL_II_HOST_DEVICE void
      integrate_value(ViewType u);

      /**
       * Evaluate the gradients of the finite element function at the quadrature
       * points.
       */
      template <typename ViewTypeIn, typename ViewTypeOut>
      DEAL_II_HOST_DEVICE void
      gradient_at_quad_pts(const ViewTypeIn u, ViewTypeOut grad_u);

      /**
       * Evaluate the values and the gradients of the finite element function at
       * the quadrature points.
       */
      template <typename ViewType1, typename ViewType2>
      DEAL_II_HOST_DEVICE void
      value_and_gradient_at_quad_pts(ViewType1 u, ViewType2 grad_u);

      /**
       * Helper function for integrate(). Integrate the gradients of the finite
       * element function.
       */
      template <bool add, typename ViewType1, typename ViewType2>
      DEAL_II_HOST_DEVICE void
      integrate_gradient(ViewType1 u, ViewType2 grad_u);

      /**
       * Helper function for integrate(). Integrate the values and the gradients
       * of the finite element function.
       */
      template <typename ViewType1, typename ViewType2>
      DEAL_II_HOST_DEVICE void
      integrate_value_and_gradient(ViewType1 u, ViewType2 grad_u);

      /**
       * TeamPolicy handle.
       */
      const TeamHandle &team_member;

      /**
       * Values of the shape functions.
       */
      Kokkos::View<Number *, MemorySpace::Default::kokkos_space> shape_values;

      /**
       * Values of the shape function gradients.
       */
      Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
        shape_gradients;

      /**
       * Values of the shape function gradients for collocation methods.
       */
      Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
        co_shape_gradients;
    };



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    DEAL_II_HOST_DEVICE
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::
      EvaluatorTensorProduct(
        const TeamHandle                                          &team_member,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space> shape_values,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
          shape_gradients,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
          co_shape_gradients)
      : team_member(team_member)
      , shape_values(shape_values)
      , shape_gradients(shape_gradients)
      , co_shape_gradients(co_shape_gradients)
    {}



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    template <int  direction,
              bool dof_to_quad,
              bool add,
              bool in_place,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::values(const ViewTypeIn in,
                                           ViewTypeOut      out) const
    {
      apply<dim, n_q_points_1d, Number, direction, dof_to_quad, add, in_place>(
        team_member, shape_values, in, out);
    }



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    template <int  direction,
              bool dof_to_quad,
              bool add,
              bool in_place,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::gradients(const ViewTypeIn in,
                                              ViewTypeOut      out) const
    {
      apply<dim, n_q_points_1d, Number, direction, dof_to_quad, add, in_place>(
        team_member, shape_gradients, in, out);
    }



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    template <int  direction,
              bool dof_to_quad,
              bool add,
              bool in_place,
              typename ViewTypeIn,
              typename ViewTypeOut>
    DEAL_II_HOST_DEVICE void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::co_gradients(const ViewTypeIn in,
                                                 ViewTypeOut      out) const
    {
      apply<dim, n_q_points_1d, Number, direction, dof_to_quad, add, in_place>(
        team_member, co_shape_gradients, in, out);
    }



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    template <typename ViewType>
    DEAL_II_HOST_DEVICE inline void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::value_at_quad_pts(ViewType u)
    {
      switch (dim)
        {
          case 1:
            {
              values<0, true, false, true>(u, u);

              break;
            }
          case 2:
            {
              values<0, true, false, true>(u, u);
              team_member.team_barrier();
              values<1, true, false, true>(u, u);

              break;
            }
          case 3:
            {
              values<0, true, false, true>(u, u);
              team_member.team_barrier();
              values<1, true, false, true>(u, u);
              team_member.team_barrier();
              values<2, true, false, true>(u, u);

              break;
            }
          default:
            {
              // Do nothing. We should throw but we can't from a __device__
              // function.
            }
        }
    }



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    template <typename ViewType>
    DEAL_II_HOST_DEVICE inline void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::integrate_value(ViewType u)
    {
      switch (dim)
        {
          case 1:
            {
              values<0, false, false, true>(u, u);

              break;
            }
          case 2:
            {
              values<0, false, false, true>(u, u);
              team_member.team_barrier();
              values<1, false, false, true>(u, u);

              break;
            }
          case 3:
            {
              values<0, false, false, true>(u, u);
              team_member.team_barrier();
              values<1, false, false, true>(u, u);
              team_member.team_barrier();
              values<2, false, false, true>(u, u);

              break;
            }
          default:
            {
              // Do nothing. We should throw but we can't from a __device__
              // function.
            }
        }
    }



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    template <typename ViewTypeIn, typename ViewTypeOut>
    DEAL_II_HOST_DEVICE inline void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::gradient_at_quad_pts(const ViewTypeIn u,
                                                         ViewTypeOut grad_u)
    {
      switch (dim)
        {
          case 1:
            {
              gradients<0, true, false, false>(
                u, Kokkos::subview(grad_u, Kokkos::ALL, 0));

              break;
            }
          case 2:
            {
              gradients<0, true, false, false>(
                u, Kokkos::subview(grad_u, Kokkos::ALL, 0));
              values<0, true, false, false>(
                u, Kokkos::subview(grad_u, Kokkos::ALL, 1));

              team_member.team_barrier();

              values<1, true, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 0),
                Kokkos::subview(grad_u, Kokkos::ALL, 0));
              gradients<1, true, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 1),
                Kokkos::subview(grad_u, Kokkos::ALL, 1));

              break;
            }
          case 3:
            {
              gradients<0, true, false, false>(
                u, Kokkos::subview(grad_u, Kokkos::ALL, 0));
              values<0, true, false, false>(
                u, Kokkos::subview(grad_u, Kokkos::ALL, 1));
              values<0, true, false, false>(
                u, Kokkos::subview(grad_u, Kokkos::ALL, 2));

              team_member.team_barrier();

              values<1, true, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 0),
                Kokkos::subview(grad_u, Kokkos::ALL, 0));
              gradients<1, true, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 1),
                Kokkos::subview(grad_u, Kokkos::ALL, 1));
              values<1, true, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 2),
                Kokkos::subview(grad_u, Kokkos::ALL, 2));

              team_member.team_barrier();

              values<2, true, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 0),
                Kokkos::subview(grad_u, Kokkos::ALL, 0));
              values<2, true, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 1),
                Kokkos::subview(grad_u, Kokkos::ALL, 1));
              gradients<2, true, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 2),
                Kokkos::subview(grad_u, Kokkos::ALL, 2));

              break;
            }
          default:
            {
              // Do nothing. We should throw but we can't from a __device__
              // function.
            }
        }
    }



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    template <typename ViewType1, typename ViewType2>
    DEAL_II_HOST_DEVICE inline void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::value_and_gradient_at_quad_pts(ViewType1 u,
                                                                   ViewType2
                                                                     grad_u)
    {
      switch (dim)
        {
          case 1:
            {
              values<0, true, false, true>(u, u);
              team_member.team_barrier();

              co_gradients<0, true, false, false>(
                u, Kokkos::subview(grad_u, Kokkos::ALL, 0));

              break;
            }
          case 2:
            {
              values<0, true, false, true>(u, u);
              team_member.team_barrier();
              values<1, true, false, true>(u, u);
              team_member.team_barrier();

              co_gradients<0, true, false, false>(
                u, Kokkos::subview(grad_u, Kokkos::ALL, 0));
              co_gradients<1, true, false, false>(
                u, Kokkos::subview(grad_u, Kokkos::ALL, 1));

              break;
            }
          case 3:
            {
              values<0, true, false, true>(u, u);
              team_member.team_barrier();
              values<1, true, false, true>(u, u);
              team_member.team_barrier();
              values<2, true, false, true>(u, u);
              team_member.team_barrier();

              co_gradients<0, true, false, false>(
                u, Kokkos::subview(grad_u, Kokkos::ALL, 0));
              co_gradients<1, true, false, false>(
                u, Kokkos::subview(grad_u, Kokkos::ALL, 1));
              co_gradients<2, true, false, false>(
                u, Kokkos::subview(grad_u, Kokkos::ALL, 2));

              break;
            }
          default:
            {
              // Do nothing. We should throw but we can't from a __device__
              // function.
            }
        }
    }



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    template <bool add, typename ViewType1, typename ViewType2>
    DEAL_II_HOST_DEVICE inline void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::integrate_gradient(ViewType1 u,
                                                       ViewType2 grad_u)
    {
      switch (dim)
        {
          case 1:
            {
              gradients<0, false, add, false>(
                Kokkos::subview(grad_u, Kokkos::ALL, dim), u);

              break;
            }
          case 2:
            {
              gradients<0, false, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 0),
                Kokkos::subview(grad_u, Kokkos::ALL, 0));
              values<0, false, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 1),
                Kokkos::subview(grad_u, Kokkos::ALL, 1));

              team_member.team_barrier();

              values<1, false, add, false>(
                Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
              team_member.team_barrier();
              gradients<1, false, true, false>(
                Kokkos::subview(grad_u, Kokkos::ALL, 1), u);

              break;
            }
          case 3:
            {
              gradients<0, false, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 0),
                Kokkos::subview(grad_u, Kokkos::ALL, 0));
              values<0, false, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 1),
                Kokkos::subview(grad_u, Kokkos::ALL, 1));
              values<0, false, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 2),
                Kokkos::subview(grad_u, Kokkos::ALL, 2));

              team_member.team_barrier();

              values<1, false, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 0),
                Kokkos::subview(grad_u, Kokkos::ALL, 0));
              gradients<1, false, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 1),
                Kokkos::subview(grad_u, Kokkos::ALL, 1));
              values<1, false, false, true>(
                Kokkos::subview(grad_u, Kokkos::ALL, 2),
                Kokkos::subview(grad_u, Kokkos::ALL, 2));

              team_member.team_barrier();

              values<2, false, add, false>(
                Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
              team_member.team_barrier();
              values<2, false, true, false>(
                Kokkos::subview(grad_u, Kokkos::ALL, 1), u);
              team_member.team_barrier();
              gradients<2, false, true, false>(
                Kokkos::subview(grad_u, Kokkos::ALL, 2), u);

              break;
            }
          default:
            {
              // Do nothing. We should throw but we can't from a __device__
              // function.
            }
        }
    }



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    template <typename ViewType1, typename ViewType2>
    DEAL_II_HOST_DEVICE inline void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::integrate_value_and_gradient(ViewType1 u,
                                                                 ViewType2
                                                                   grad_u)
    {
      switch (dim)
        {
          case 1:
            {
              co_gradients<0, false, true, false>(
                Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
              team_member.team_barrier();

              values<0, false, false, true>(u, u);

              break;
            }
          case 2:
            {
              co_gradients<1, false, true, false>(
                Kokkos::subview(grad_u, Kokkos::ALL, 1), u);
              team_member.team_barrier();
              co_gradients<0, false, true, false>(
                Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
              team_member.team_barrier();

              values<1, false, false, true>(u, u);
              team_member.team_barrier();
              values<0, false, false, true>(u, u);
              team_member.team_barrier();

              break;
            }
          case 3:
            {
              co_gradients<2, false, true, false>(
                Kokkos::subview(grad_u, Kokkos::ALL, 2), u);
              team_member.team_barrier();
              co_gradients<1, false, true, false>(
                Kokkos::subview(grad_u, Kokkos::ALL, 1), u);
              team_member.team_barrier();
              co_gradients<0, false, true, false>(
                Kokkos::subview(grad_u, Kokkos::ALL, 0), u);
              team_member.team_barrier();

              values<2, false, false, true>(u, u);
              team_member.team_barrier();
              values<1, false, false, true>(u, u);
              team_member.team_barrier();
              values<0, false, false, true>(u, u);
              team_member.team_barrier();

              break;
            }
          default:
            {
              // Do nothing. We should throw but we can't from a __device__
              // function.
            }
        }
    }
  } // namespace internal
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
