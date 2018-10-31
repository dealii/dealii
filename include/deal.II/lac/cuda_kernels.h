// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

#ifndef dealii_cuda_kernels_h
#define dealii_cuda_kernels_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_COMPILER_CUDA_AWARE


#  include <deal.II/base/cuda_size.h>
#  include <deal.II/base/types.h>

#  include <deal.II/lac/cuda_atomic.h>

#  include <assert.h>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace CUDAWrappers
  {
    /**
     * Namespace containing the CUDA kernels.
     */
    namespace kernel
    {
      using ::dealii::CUDAWrappers::block_size;
      using ::dealii::CUDAWrappers::chunk_size;
      using size_type = types::global_dof_index;

      /**
       * Multiply each entry of @p val of size @p N by @p a.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      __global__ void
      vec_scale(Number *val, const Number a, const size_type N);



      /**
       * Functor defining the addition of two Numbers.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      struct Binop_Addition
      {
        __device__ static inline Number
        operation(const Number a, const Number b)
        {
          return a + b;
        }
      };

      template <typename Number>
      struct Binop_Addition<std::complex<Number>>
      {
        __device__ static inline std::complex<Number>
        operation(const std::complex<Number> a, const std::complex<Number> b)
        {
          printf("This function is not implemented for std::complex<Number>!");
          assert(false);
          return {};
        }
      };



      /**
       * Functor defining the subtraction of two Numbers.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      struct Binop_Subtraction
      {
        __device__ static inline Number
        operation(const Number a, const Number b)
        {
          return a - b;
        }
      };



      /**
       * Apply the functor @tparam Binop to each element of @p v1 and @p v2.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number, template <typename> class Binop>
      __global__ void
      vector_bin_op(Number *v1, const Number *v2, const size_type N);



      /**
       * Structure implementing the functions used to add elements when using a
       * reduction.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      struct ElemSum
      {
        __device__ static Number
        reduction_op(const Number a, const Number b);

        __device__ static Number
        atomic_op(Number *dst, const Number a);

        __device__ static Number
        element_wise_op(const Number a);

        __device__ static Number
        null_value();
      };



      /**
       * Structure implementing the functions used to compute the L1 norm when
       * using a reduction.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      struct L1Norm
      {
        __device__ static Number
        reduction_op(const Number a, const Number b);

        __device__ static Number
        atomic_op(Number *dst, const Number a);

        __device__ static Number
        element_wise_op(const Number a);

        __device__ static Number
        null_value();
      };



      /**
       * Structure implementing the functions used to compute the L-infinity
       * norm when using a reduction.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      struct LInfty
      {
        __device__ static Number
        reduction_op(const Number a, const Number b);

        __device__ static Number
        atomic_op(Number *dst, const Number a);

        __device__ static Number
        element_wise_op(const Number a);

        __device__ static Number
        null_value();
      };



      /**
       * Perform a reduction on @p v using @tparam Operation.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number, typename Operation>
      __global__ void
      reduction(Number *result, const Number *v, const size_type N);



      /**
       * Structure implementing the functions used to compute the dot product
       * norm when using a double vector reduction.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      struct DotProduct
      {
        __device__ static Number
        binary_op(const Number a, const Number b);

        __device__ static Number
        reduction_op(const Number a, const Number b);

        __device__ static Number
        atomic_op(Number *dst, const Number a);

        __device__ static Number
        null_value();
      };



      /**
       * Perform a binary operation on each element of @p v1 and @p v2 followed
       * by reduction on the resulting array.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number, typename Operation>
      __global__ void
      double_vector_reduction(Number *        result,
                              const Number *  v1,
                              const Number *  v2,
                              const size_type N);



      /**
       * Add @p a to each element of @p val.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      __global__ void
      vec_add(Number *val, const Number a, const size_type N);



      /**
       * Addition of a multiple of a vector, i.e., <tt>val += a*V_val</tt>.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      __global__ void
      add_aV(Number *        val,
             const Number    a,
             const Number *  V_val,
             const size_type N);



      /**
       * Addition of multiple scaled vector, i.e., <tt>val += a*V_val +
       * b*W_val</tt>.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      __global__ void
      add_aVbW(Number *        val,
               const Number    a,
               const Number *  V_val,
               const Number    b,
               const Number *  W_val,
               const size_type N);



      /**
       * Scaling and simple addition of a multiple of a vector, i.e. <tt>val =
       * = s*val + a*V_val</tt>
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      __global__ void
      sadd(const Number    s,
           Number *        val,
           const Number    a,
           const Number *  V_val,
           const size_type N);



      /**
       * Scaling and multiple additions of scaled vectors, i.e. <tt>val =
       * = s*val + a*V_val + b*W_val</tt>
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      __global__ void
      sadd(const Number    s,
           Number *        val,
           const Number    a,
           const Number *  V_val,
           const Number    b,
           const Number *  W_val,
           const size_type N);



      /**
       * Scale each element of this vector by the corresponding element in the
       * argument.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      __global__ void
      scale(Number *val, const Number *V_val, const size_type N);



      /**
       * Assignment <tt>val = a*V_val</tt>.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      __global__ void
      equ(Number *val, const Number a, const Number *V_val, const size_type N);



      /**
       * Assignment <tt>val = a*V_val + b*W_val</tt>.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      __global__ void
      equ(Number *        val,
          const Number    a,
          const Number *  V_val,
          const Number    b,
          const Number *  W_val,
          const size_type N);



      /**
       * Perform a combined operation of a vector addition and a subsequent
       * inner product, returning the value of the inner product.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      __global__ void
      add_and_dot(Number *        res,
                  Number *        v1,
                  const Number *  v2,
                  const Number *  v3,
                  const Number    a,
                  const size_type N);



      /**
       * Set each element of @p val to @p s.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      __global__ void
      set(Number *val, const Number s, const size_type N);


      /**
       * Set each element @v val to @p v using @p indices as permutation, i.e.,
       * <tt>val[indices[i]] = v[i]</tt>.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      __global__ void
      set_permutated(Number *         val,
                     const Number *   v,
                     const size_type *indices,
                     const size_type  N);



      /**
       * Add each element @v val to @p v using @p indices as permutation, i.e.,
       * <tt>val[indices[i]] += v[i]</tt>.
       *
       * @ingroup CUDAWrappers
       */
      template <typename Number>
      __global__ void
      add_permutated(Number *         val,
                     const Number *   v,
                     const size_type *indices,
                     const size_type  N);
    } // namespace kernel
  }   // namespace CUDAWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
