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


#ifndef dealii_cuda_tensor_product_kernels_h
#define dealii_cuda_tensor_product_kernels_h

#include <deal.II/base/config.h>

#include <deal.II/base/utilities.h>

#include <deal.II/matrix_free/cuda_matrix_free.templates.h>

#ifdef DEAL_II_WITH_CUDA

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

    /**
     * Helper function for values() and gradients().
     */
    template <int dim,
              int n_q_points_1d,
              typename Number,
              int  direction,
              bool dof_to_quad,
              bool add,
              bool in_place>
    DEAL_II_HOST_DEVICE void
    apply(Kokkos::View<Number *, MemorySpace::Default::kokkos_space> shape_data,
          const Number *                                             in,
          Number *                                                   out)
    {
      KOKKOS_IF_ON_DEVICE(
        const unsigned int i = (dim == 1) ? 0 : threadIdx.x % n_q_points_1d;
        const unsigned int j = (dim == 3) ? threadIdx.y : 0;
        const unsigned int q = (dim == 1) ? (threadIdx.x % n_q_points_1d) :
                               (dim == 2) ? threadIdx.y :
                                            threadIdx.z;

        // This loop simply multiply the shape function at the quadrature point
        // by the value finite element coefficient.
        Number t = 0;
        for (int k = 0; k < n_q_points_1d; ++k) {
          const unsigned int shape_idx =
            dof_to_quad ? (q + k * n_q_points_1d) : (k + q * n_q_points_1d);
          const unsigned int source_idx =
            (direction == 0) ? (k + n_q_points_1d * (i + n_q_points_1d * j)) :
            (direction == 1) ? (i + n_q_points_1d * (k + n_q_points_1d * j)) :
                               (i + n_q_points_1d * (j + n_q_points_1d * k));
          t += shape_data[shape_idx] *
               (in_place ? out[source_idx] : in[source_idx]);
        }

        if (in_place) __syncthreads();

        const unsigned int destination_idx =
          (direction == 0) ? (q + n_q_points_1d * (i + n_q_points_1d * j)) :
          (direction == 1) ? (i + n_q_points_1d * (q + n_q_points_1d * j)) :
                             (i + n_q_points_1d * (j + n_q_points_1d * q));

        if (add) out[destination_idx] += t;
        else out[destination_idx] = t;)
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
      DEAL_II_HOST_DEVICE
      EvaluatorTensorProduct(
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space> shape_values,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
          shape_gradients,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
          co_shape_gradients);

      /**
       * Evaluate the values of a finite element function at the quadrature
       * points.
       */
      template <int direction, bool dof_to_quad, bool add, bool in_place>
      DEAL_II_HOST_DEVICE void
      values(const Number *in, Number *out) const;

      /**
       * Evaluate the gradient of a finite element function at the quadrature
       * points for a given @p direction.
       */
      template <int direction, bool dof_to_quad, bool add, bool in_place>
      DEAL_II_HOST_DEVICE void
      gradients(const Number *in, Number *out) const;

      /**
       * TODO
       */
      template <int direction, bool dof_to_quad, bool add, bool in_place>
      DEAL_II_HOST_DEVICE void
      co_gradients(const Number *in, Number *out) const;

      /**
       * Evaluate the finite element function at the quadrature points.
       */
      DEAL_II_HOST_DEVICE void
      value_at_quad_pts(Number *u);

      /**
       * Helper function for integrate(). Integrate the finite element function.
       */
      DEAL_II_HOST_DEVICE void
      integrate_value(Number *u);

      /**
       * Evaluate the gradients of the finite element function at the quadrature
       * points.
       */
      DEAL_II_HOST_DEVICE void
      gradient_at_quad_pts(const Number *const u, Number *grad_u[dim]);

      /**
       * Evaluate the values and the gradients of the finite element function at
       * the quadrature points.
       */
      DEAL_II_HOST_DEVICE void
      value_and_gradient_at_quad_pts(Number *const u, Number *grad_u[dim]);

      /**
       * Helper function for integrate(). Integrate the gradients of the finite
       * element function.
       */
      template <bool add>
      DEAL_II_HOST_DEVICE void
      integrate_gradient(Number *u, Number *grad_u[dim]);

      /**
       * Helper function for integrate(). Integrate the values and the gradients
       * of the finite element function.
       */
      DEAL_II_HOST_DEVICE void
      integrate_value_and_gradient(Number *u, Number *grad_u[dim]);

      // TODO shape values
      Kokkos::View<Number *, MemorySpace::Default::kokkos_space> shape_values;

      // TODO shape gradients
      Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
        shape_gradients;

      // TODO shape gradients for collocation methods
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
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space> shape_values,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
          shape_gradients,
        Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
          co_shape_gradients)
      : shape_values(shape_values)
      , shape_gradients(shape_gradients)
      , co_shape_gradients(co_shape_gradients)
    {}



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    template <int direction, bool dof_to_quad, bool add, bool in_place>
    DEAL_II_HOST_DEVICE void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::values(const Number *in, Number *out) const
    {
      apply<dim, n_q_points_1d, Number, direction, dof_to_quad, add, in_place>(
        shape_values, in, out);
    }



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    template <int direction, bool dof_to_quad, bool add, bool in_place>
    DEAL_II_HOST_DEVICE void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::gradients(const Number *in,
                                              Number *      out) const
    {
      apply<dim, n_q_points_1d, Number, direction, dof_to_quad, add, in_place>(
        shape_gradients, in, out);
    }



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    template <int direction, bool dof_to_quad, bool add, bool in_place>
    DEAL_II_HOST_DEVICE void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::co_gradients(const Number *in,
                                                 Number *      out) const
    {
      apply<dim, n_q_points_1d, Number, direction, dof_to_quad, add, in_place>(
        co_shape_gradients, in, out);
    }



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    DEAL_II_HOST_DEVICE inline void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::value_at_quad_pts(Number *u)
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
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              values<1, true, false, true>(u, u);

              break;
            }
          case 3:
            {
              values<0, true, false, true>(u, u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              values<1, true, false, true>(u, u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
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
    DEAL_II_HOST_DEVICE inline void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::integrate_value(Number *u)
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
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              values<1, false, false, true>(u, u);

              break;
            }
          case 3:
            {
              values<0, false, false, true>(u, u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              values<1, false, false, true>(u, u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
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
    DEAL_II_HOST_DEVICE inline void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::gradient_at_quad_pts(const Number *const u,
                                                         Number *grad_u[dim])
    {
      switch (dim)
        {
          case 1:
            {
              gradients<0, true, false, false>(u, grad_u[0]);

              break;
            }
          case 2:
            {
              gradients<0, true, false, false>(u, grad_u[0]);
              values<0, true, false, false>(u, grad_u[1]);

              KOKKOS_IF_ON_DEVICE(__syncthreads();)

              values<1, true, false, true>(grad_u[0], grad_u[0]);
              gradients<1, true, false, true>(grad_u[1], grad_u[1]);

              break;
            }
          case 3:
            {
              gradients<0, true, false, false>(u, grad_u[0]);
              values<0, true, false, false>(u, grad_u[1]);
              values<0, true, false, false>(u, grad_u[2]);

              KOKKOS_IF_ON_DEVICE(__syncthreads();)

              values<1, true, false, true>(grad_u[0], grad_u[0]);
              gradients<1, true, false, true>(grad_u[1], grad_u[1]);
              values<1, true, false, true>(grad_u[2], grad_u[2]);

              KOKKOS_IF_ON_DEVICE(__syncthreads();)

              values<2, true, false, true>(grad_u[0], grad_u[0]);
              values<2, true, false, true>(grad_u[1], grad_u[1]);
              gradients<2, true, false, true>(grad_u[2], grad_u[2]);

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
    DEAL_II_HOST_DEVICE inline void
    EvaluatorTensorProduct<
      evaluate_general,
      dim,
      fe_degree,
      n_q_points_1d,
      Number>::value_and_gradient_at_quad_pts(Number *const u,
                                              Number *      grad_u[dim])
    {
      switch (dim)
        {
          case 1:
            {
              values<0, true, false, true>(u, u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)

              co_gradients<0, true, false, false>(u, grad_u[0]);

              break;
            }
          case 2:
            {
              values<0, true, false, true>(u, u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              values<1, true, false, true>(u, u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)

              co_gradients<0, true, false, false>(u, grad_u[0]);
              co_gradients<1, true, false, false>(u, grad_u[1]);

              break;
            }
          case 3:
            {
              values<0, true, false, true>(u, u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              values<1, true, false, true>(u, u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              values<2, true, false, true>(u, u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)

              co_gradients<0, true, false, false>(u, grad_u[0]);
              co_gradients<1, true, false, false>(u, grad_u[1]);
              co_gradients<2, true, false, false>(u, grad_u[2]);

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
    template <bool add>
    DEAL_II_HOST_DEVICE inline void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::integrate_gradient(Number *u,
                                                       Number *grad_u[dim])
    {
      switch (dim)
        {
          case 1:
            {
              gradients<0, false, add, false>(

                grad_u[dim], u);

              break;
            }
          case 2:
            {
              gradients<0, false, false, true>(grad_u[0], grad_u[0]);
              values<0, false, false, true>(grad_u[1], grad_u[1]);

              KOKKOS_IF_ON_DEVICE(__syncthreads();)

              values<1, false, add, false>(grad_u[0], u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              gradients<1, false, true, false>(grad_u[1], u);

              break;
            }
          case 3:
            {
              gradients<0, false, false, true>(grad_u[0], grad_u[0]);
              values<0, false, false, true>(grad_u[1], grad_u[1]);
              values<0, false, false, true>(grad_u[2], grad_u[2]);

              KOKKOS_IF_ON_DEVICE(__syncthreads();)

              values<1, false, false, true>(grad_u[0], grad_u[0]);
              gradients<1, false, false, true>(grad_u[1], grad_u[1]);
              values<1, false, false, true>(grad_u[2], grad_u[2]);

              KOKKOS_IF_ON_DEVICE(__syncthreads();)

              values<2, false, add, false>(grad_u[0], u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              values<2, false, true, false>(grad_u[1], u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              gradients<2, false, true, false>(grad_u[2], u);

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
    DEAL_II_HOST_DEVICE inline void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::integrate_value_and_gradient(Number *u,
                                                                 Number
                                                                   *grad_u[dim])
    {
      switch (dim)
        {
          case 1:
            {
              co_gradients<0, false, true, false>(grad_u[0], u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)

              values<0, false, false, true>(u, u);

              break;
            }
          case 2:
            {
              co_gradients<1, false, true, false>(grad_u[1], u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              co_gradients<0, false, true, false>(grad_u[0], u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)

              values<1, false, false, true>(u, u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              values<0, false, false, true>(u, u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)

              break;
            }
          case 3:
            {
              co_gradients<2, false, true, false>(grad_u[2], u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              co_gradients<1, false, true, false>(grad_u[1], u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              co_gradients<0, false, true, false>(grad_u[0], u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)

              values<2, false, false, true>(u, u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              values<1, false, false, true>(u, u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)
              values<0, false, false, true>(u, u);
              KOKKOS_IF_ON_DEVICE(__syncthreads();)

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

#endif
