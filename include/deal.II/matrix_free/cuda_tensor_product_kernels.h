// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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

#ifdef DEAL_II_COMPILER_CUDA_AWARE

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
      static constexpr unsigned int dofs_per_cell =
        Utilities::pow(fe_degree + 1, dim);
      static constexpr unsigned int n_q_points =
        Utilities::pow(n_q_points_1d, dim);

      __device__
      EvaluatorTensorProduct();

      /**
       * Evaluate the values of a finite element function at the quadrature
       * points.
       */
      template <int direction, bool dof_to_quad, bool add, bool in_place>
      __device__ void
      values(Number shape_values[], const Number *in, Number *out) const;

      /**
       * Evaluate the gradient of a finite element function at the quadrature
       * points for a given @p direction.
       */
      template <int direction, bool dof_to_quad, bool add, bool in_place>
      __device__ void
      gradients(Number shape_gradients[], const Number *in, Number *out) const;

      /**
       * Helper function for values() and gradients().
       */
      template <int direction, bool dof_to_quad, bool add, bool in_place>
      __device__ void
      apply(Number shape_data[], const Number *in, Number *out) const;

      /**
       * Evaluate the finite element function at the quadrature points.
       */
      __device__ void
      value_at_quad_pts(Number *u);

      /**
       * Helper function for integrate(). Integrate the finite element function.
       */
      __device__ void
      integrate_value(Number *u);

      /**
       * Evaluate the gradients of the finite element function at the quadrature
       * points.
       */
      __device__ void
      gradient_at_quad_pts(const Number *const u, Number *grad_u[dim]);

      /**
       * Evaluate the values and the gradients of the finite element function at
       *  the quadrature points.
       */
      __device__ void
      value_and_gradient_at_quad_pts(Number *const u, Number *grad_u[dim]);

      /**
       * Helper function for integrate(). Integrate the gradients of the finite
       * element function.
       */
      template <bool add>
      __device__ void
      integrate_gradient(Number *u, Number *grad_u[dim]);

      /**
       * Helper function for integrate(). Integrate the values and the gradients
       * of the finite element function.
       */
      __device__ void
      integrate_value_and_gradient(Number *u, Number *grad_u[dim]);
    };



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    __device__
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::EvaluatorTensorProduct()
    {}



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    template <int direction, bool dof_to_quad, bool add, bool in_place>
    __device__ void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::values(Number        shape_values[],
                                           const Number *in,
                                           Number *      out) const
    {
      apply<direction, dof_to_quad, add, in_place>(shape_values, in, out);
    }



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    template <int direction, bool dof_to_quad, bool add, bool in_place>
    __device__ void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::gradients(Number        shape_gradients[],
                                              const Number *in,
                                              Number *      out) const
    {
      apply<direction, dof_to_quad, add, in_place>(shape_gradients, in, out);
    }



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    template <int direction, bool dof_to_quad, bool add, bool in_place>
    __device__ void
    EvaluatorTensorProduct<evaluate_general,
                           dim,
                           fe_degree,
                           n_q_points_1d,
                           Number>::apply(Number        shape_data[],
                                          const Number *in,
                                          Number *      out) const
    {
      const unsigned int i = (dim == 1) ? 0 : threadIdx.x % n_q_points_1d;
      const unsigned int j = (dim == 3) ? threadIdx.y : 0;
      const unsigned int q = (dim == 1) ?
                               (threadIdx.x % n_q_points_1d) :
                               (dim == 2) ? threadIdx.y : threadIdx.z;

      // This loop simply multiply the shape function at the quadrature point by
      // the value finite element coefficient.
      Number t = 0;
      for (int k = 0; k < n_q_points_1d; ++k)
        {
          const unsigned int shape_idx =
            dof_to_quad ? (q + k * n_q_points_1d) : (k + q * n_q_points_1d);
          const unsigned int source_idx =
            (direction == 0) ?
              (k + n_q_points_1d * (i + n_q_points_1d * j)) :
              (direction == 1) ? (i + n_q_points_1d * (k + n_q_points_1d * j)) :
                                 (i + n_q_points_1d * (j + n_q_points_1d * k));
          t += shape_data[shape_idx] *
               (in_place ? out[source_idx] : in[source_idx]);
        }

      if (in_place)
        __syncthreads();

      const unsigned int destination_idx =
        (direction == 0) ?
          (q + n_q_points_1d * (i + n_q_points_1d * j)) :
          (direction == 1) ? (i + n_q_points_1d * (q + n_q_points_1d * j)) :
                             (i + n_q_points_1d * (j + n_q_points_1d * q));

      if (add)
        out[destination_idx] += t;
      else
        out[destination_idx] = t;
    }



    template <int dim, int fe_degree, int n_q_points_1d, typename Number>
    inline __device__ void
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
              values<0, true, false, true>(global_shape_values, u, u);

              break;
            }
          case 2:
            {
              values<0, true, false, true>(global_shape_values, u, u);
              __syncthreads();
              values<1, true, false, true>(global_shape_values, u, u);

              break;
            }
          case 3:
            {
              values<0, true, false, true>(global_shape_values, u, u);
              __syncthreads();
              values<1, true, false, true>(global_shape_values, u, u);
              __syncthreads();
              values<2, true, false, true>(global_shape_values, u, u);

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
    inline __device__ void
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
              values<0, false, false, true>(global_shape_values, u, u);

              break;
            }
          case 2:
            {
              values<0, false, false, true>(global_shape_values, u, u);
              __syncthreads();
              values<1, false, false, true>(global_shape_values, u, u);

              break;
            }
          case 3:
            {
              values<0, false, false, true>(global_shape_values, u, u);
              __syncthreads();
              values<1, false, false, true>(global_shape_values, u, u);
              __syncthreads();
              values<2, false, false, true>(global_shape_values, u, u);

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
    inline __device__ void
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
              gradients<0, true, false, false>(global_shape_gradients,
                                               u,
                                               grad_u[0]);

              break;
            }
          case 2:
            {
              gradients<0, true, false, false>(global_shape_gradients,
                                               u,
                                               grad_u[0]);
              values<0, true, false, false>(global_shape_values, u, grad_u[1]);

              __syncthreads();

              values<1, true, false, true>(global_shape_values,
                                           grad_u[0],
                                           grad_u[0]);
              gradients<1, true, false, true>(global_shape_gradients,
                                              grad_u[1],
                                              grad_u[1]);

              break;
            }
          case 3:
            {
              gradients<0, true, false, false>(global_shape_gradients,
                                               u,
                                               grad_u[0]);
              values<0, true, false, false>(global_shape_values, u, grad_u[1]);
              values<0, true, false, false>(global_shape_values, u, grad_u[2]);

              __syncthreads();

              values<1, true, false, true>(global_shape_values,
                                           grad_u[0],
                                           grad_u[0]);
              gradients<1, true, false, true>(global_shape_gradients,
                                              grad_u[1],
                                              grad_u[1]);
              values<1, true, false, true>(global_shape_values,
                                           grad_u[2],
                                           grad_u[2]);

              __syncthreads();

              values<2, true, false, true>(global_shape_values,
                                           grad_u[0],
                                           grad_u[0]);
              values<2, true, false, true>(global_shape_values,
                                           grad_u[1],
                                           grad_u[1]);
              gradients<2, true, false, true>(global_shape_gradients,
                                              grad_u[2],
                                              grad_u[2]);

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
    inline __device__ void
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
              values<0, true, false, true>(global_shape_values, u, u);
              __syncthreads();

              gradients<0, true, false, false>(global_co_shape_gradients,
                                               u,
                                               grad_u[0]);

              break;
            }
          case 2:
            {
              values<0, true, false, true>(global_shape_values, u, u);
              __syncthreads();
              values<1, true, false, true>(global_shape_values, u, u);
              __syncthreads();

              gradients<0, true, false, false>(global_co_shape_gradients,
                                               u,
                                               grad_u[0]);
              gradients<1, true, false, false>(global_co_shape_gradients,
                                               u,
                                               grad_u[1]);

              break;
            }
          case 3:
            {
              values<0, true, false, true>(global_shape_values, u, u);
              __syncthreads();
              values<1, true, false, true>(global_shape_values, u, u);
              __syncthreads();
              values<2, true, false, true>(global_shape_values, u, u);
              __syncthreads();

              gradients<0, true, false, false>(global_co_shape_gradients,
                                               u,
                                               grad_u[0]);
              gradients<1, true, false, false>(global_co_shape_gradients,
                                               u,
                                               grad_u[1]);
              gradients<2, true, false, false>(global_co_shape_gradients,
                                               u,
                                               grad_u[2]);

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
    inline __device__ void
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
              gradients<0, false, add, false>(global_shape_gradients,
                                              grad_u[dim],
                                              u);

              break;
            }
          case 2:
            {
              gradients<0, false, false, true>(global_shape_gradients,
                                               grad_u[0],
                                               grad_u[0]);
              values<0, false, false, true>(global_shape_values,
                                            grad_u[1],
                                            grad_u[1]);

              __syncthreads();

              values<1, false, add, false>(global_shape_values, grad_u[0], u);
              __syncthreads();
              gradients<1, false, true, false>(global_shape_gradients,
                                               grad_u[1],
                                               u);

              break;
            }
          case 3:
            {
              gradients<0, false, false, true>(global_shape_gradients,
                                               grad_u[0],
                                               grad_u[0]);
              values<0, false, false, true>(global_shape_values,
                                            grad_u[1],
                                            grad_u[1]);
              values<0, false, false, true>(global_shape_values,
                                            grad_u[2],
                                            grad_u[2]);

              __syncthreads();

              values<1, false, false, true>(global_shape_values,
                                            grad_u[0],
                                            grad_u[0]);
              gradients<1, false, false, true>(global_shape_gradients,
                                               grad_u[1],
                                               grad_u[1]);
              values<1, false, false, true>(global_shape_values,
                                            grad_u[2],
                                            grad_u[2]);

              __syncthreads();

              values<2, false, add, false>(global_shape_values, grad_u[0], u);
              __syncthreads();
              values<2, false, true, false>(global_shape_values, grad_u[1], u);
              __syncthreads();
              gradients<2, false, true, false>(global_shape_gradients,
                                               grad_u[2],
                                               u);

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
    inline __device__ void
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
              gradients<0, false, true, false>(global_co_shape_gradients,
                                               grad_u[0],
                                               u);
              __syncthreads();

              values<0, false, false, true>(global_shape_values, u, u);

              break;
            }
          case 2:
            {
              gradients<1, false, true, false>(global_co_shape_gradients,
                                               grad_u[1],
                                               u);
              __syncthreads();
              gradients<0, false, true, false>(global_co_shape_gradients,
                                               grad_u[0],
                                               u);
              __syncthreads();

              values<1, false, false, true>(global_shape_values, u, u);
              __syncthreads();
              values<0, false, false, true>(global_shape_values, u, u);
              __syncthreads();

              break;
            }
          case 3:
            {
              gradients<2, false, true, false>(global_co_shape_gradients,
                                               grad_u[2],
                                               u);
              __syncthreads();
              gradients<1, false, true, false>(global_co_shape_gradients,
                                               grad_u[1],
                                               u);
              __syncthreads();
              gradients<0, false, true, false>(global_co_shape_gradients,
                                               grad_u[0],
                                               u);
              __syncthreads();

              values<2, false, false, true>(global_shape_values, u, u);
              __syncthreads();
              values<1, false, false, true>(global_shape_values, u, u);
              __syncthreads();
              values<0, false, false, true>(global_shape_values, u, u);
              __syncthreads();

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
