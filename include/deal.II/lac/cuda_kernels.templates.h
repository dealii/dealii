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

#ifndef dealii__cuda_kernels_templates_h
#define dealii__cuda_kernels_templates_h

#include <deal.II/lac/cuda_kernels.h>

#ifdef DEAL_II_COMPILER_CUDA_AWARE

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace CUDAWrappers
  {
    namespace kernel
    {
      template <typename Number>
      __global__ void
      vec_scale(Number *val, const Number a, const size_type N)
      {
        const size_type idx_base =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        for (unsigned int i = 0; i < chunk_size; ++i)
          {
            const size_type idx = idx_base + i * block_size;
            if (idx < N)
              val[idx] *= a;
          }
      }



      template <typename Number, template <typename> class Binop>
      __global__ void
      vector_bin_op(Number *v1, const Number *v2, const size_type N)
      {
        const size_type idx_base =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        for (unsigned int i = 0; i < chunk_size; ++i)
          {
            const size_type idx = idx_base + i * block_size;
            if (idx < N)
              v1[idx] = Binop<Number>::operation(v1[idx], v2[idx]);
          }
      }



      template <typename Number>
      __device__ Number
                 ElemSum<Number>::reduction_op(const Number a, const Number b)
      {
        return (a + b);
      }



      template <typename Number>
      __device__ Number
                 ElemSum<Number>::atomic_op(Number *dst, const Number a)
      {
        return atomicAdd_wrapper(dst, a);
      }



      template <typename Number>
      __device__ Number
                 ElemSum<Number>::element_wise_op(const Number a)
      {
        return a;
      }



      template <typename Number>
      __device__ Number
                 ElemSum<Number>::null_value()
      {
        return Number();
      }



      template <typename Number>
      __device__ Number
                 L1Norm<Number>::reduction_op(const Number a, const Number b)
      {
        return (a + b);
      }



      template <typename Number>
      __device__ Number
                 L1Norm<Number>::atomic_op(Number *dst, const Number a)
      {
        return atomicAdd_wrapper(dst, a);
      }



      template <typename Number>
      __device__ Number
                 L1Norm<Number>::element_wise_op(const Number a)
      {
        return std::fabs(a);
      }



      template <typename Number>
      __device__ Number
                 L1Norm<Number>::null_value()
      {
        return Number();
      }



      template <typename Number>
      __device__ Number
                 LInfty<Number>::reduction_op(const Number a, const Number b)
      {
        if (a > b)
          return a;
        else
          return b;
      }



      template <typename Number>
      __device__ Number
                 LInfty<Number>::atomic_op(Number *dst, const Number a)
      {
        return atomicMax_wrapper(dst, a);
      }



      template <typename Number>
      __device__ Number
                 LInfty<Number>::element_wise_op(const Number a)
      {
        return std::fabs(a);
      }



      template <typename Number>
      __device__ Number
                 LInfty<Number>::null_value()
      {
        return Number();
      }



      template <typename Number, typename Operation>
      __device__ void
      reduce_within_warp(volatile Number *result_buffer, size_type local_idx)
      {
        if (block_size >= 64)
          result_buffer[local_idx] =
            Operation::reduction_op(result_buffer[local_idx],
                                    result_buffer[local_idx + 32]);
        if (block_size >= 32)
          result_buffer[local_idx] =
            Operation::reduction_op(result_buffer[local_idx],
                                    result_buffer[local_idx + 16]);
        if (block_size >= 16)
          result_buffer[local_idx] =
            Operation::reduction_op(result_buffer[local_idx],
                                    result_buffer[local_idx + 8]);
        if (block_size >= 8)
          result_buffer[local_idx] =
            Operation::reduction_op(result_buffer[local_idx],
                                    result_buffer[local_idx + 4]);
        if (block_size >= 4)
          result_buffer[local_idx] =
            Operation::reduction_op(result_buffer[local_idx],
                                    result_buffer[local_idx + 2]);
        if (block_size >= 2)
          result_buffer[local_idx] =
            Operation::reduction_op(result_buffer[local_idx],
                                    result_buffer[local_idx + 1]);
      }



      template <typename Number, typename Operation>
      __device__ void
      reduce(Number *        result,
             Number *        result_buffer,
             const size_type local_idx,
             const size_type global_idx,
             const size_type N)
      {
        for (size_type s = block_size / 2; s > 32; s = s >> 1)
          {
            if (local_idx < s)
              result_buffer[local_idx] =
                Operation::reduction_op(result_buffer[local_idx],
                                        result_buffer[local_idx + s]);
            __syncthreads();
          }

        if (local_idx < 32)
          reduce_within_warp<Number, Operation>(result_buffer, local_idx);

        if (local_idx == 0)
          Operation::atomic_op(result, result_buffer[0]);
      }



      template <typename Number, typename Operation>
      __global__ void
      reduction(Number *result, const Number *v, const size_type N)
      {
        __shared__ Number result_buffer[block_size];

        const size_type global_idx =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        const size_type local_idx = threadIdx.x;

        if (global_idx < N)
          result_buffer[local_idx] = Operation::element_wise_op(v[global_idx]);
        else
          result_buffer[local_idx] = Operation::null_value();

        __syncthreads();

        reduce<Number, Operation>(
          result, result_buffer, local_idx, global_idx, N);
      }



      template <typename Number>
      __device__ Number
                 DotProduct<Number>::binary_op(const Number a, const Number b)
      {
        return a * b;
      }



      template <typename Number>
      __device__ Number
                 DotProduct<Number>::reduction_op(const Number a, const Number b)
      {
        return a + b;
      }



      template <typename Number>
      __device__ Number
                 DotProduct<Number>::atomic_op(Number *dst, const Number a)
      {
        return atomicAdd_wrapper(dst, a);
      }



      template <typename Number>
      __device__ Number
                 DotProduct<Number>::null_value()
      {
        return Number();
      }



      template <typename Number, typename Operation>
      __global__ void
      double_vector_reduction(Number *        result,
                              const Number *  v1,
                              const Number *  v2,
                              const size_type N)
      {
        __shared__ Number result_buffer[block_size];

        const size_type global_idx =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        const size_type local_idx = threadIdx.x;

        if (global_idx < N)
          result_buffer[local_idx] =
            Operation::binary_op(v1[global_idx], v2[global_idx]);
        else
          result_buffer[local_idx] = Operation::null_value();

        for (unsigned int i = 1; i < chunk_size; ++i)
          {
            const size_type idx = global_idx + i * block_size;
            if (idx < N)
              result_buffer[local_idx] =
                Operation::reduction_op(result_buffer[local_idx],
                                        Operation::binary_op(v1[idx], v2[idx]));
          }

        __syncthreads();

        reduce<Number, Operation>(
          result, result_buffer, local_idx, global_idx, N);
      }



      template <typename Number>
      __global__ void
      vec_add(Number *val, const Number a, const size_type N)
      {
        const size_type idx_base =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        for (unsigned int i = 0; i < chunk_size; ++i)
          {
            const size_type idx = idx_base + i * block_size;
            if (idx < N)
              val[idx] += a;
          }
      }



      template <typename Number>
      __global__ void
      add_aV(Number *        val,
             const Number    a,
             const Number *  V_val,
             const size_type N)
      {
        const size_type idx_base =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        for (unsigned int i = 0; i < chunk_size; ++i)
          {
            const size_type idx = idx_base + i * block_size;
            if (idx < N)
              val[idx] += a * V_val[idx];
          }
      }



      template <typename Number>
      __global__ void
      add_aVbW(Number *        val,
               const Number    a,
               const Number *  V_val,
               const Number    b,
               const Number *  W_val,
               const size_type N)
      {
        const size_type idx_base =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        for (unsigned int i = 0; i < chunk_size; ++i)
          {
            const size_type idx = idx_base + i * block_size;
            if (idx < N)
              val[idx] += a * V_val[idx] + b * W_val[idx];
          }
      }



      template <typename Number>
      __global__ void
      sadd(const Number    s,
           Number *        val,
           const Number    a,
           const Number *  V_val,
           const size_type N)
      {
        const size_type idx_base =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        for (unsigned int i = 0; i < chunk_size; ++i)
          {
            const size_type idx = idx_base + i * block_size;
            if (idx < N)
              val[idx] = s * val[idx] + a * V_val[idx];
          }
      }



      template <typename Number>
      __global__ void
      sadd(const Number    s,
           Number *        val,
           const Number    a,
           const Number *  V_val,
           const Number    b,
           const Number *  W_val,
           const size_type N)
      {
        const size_type idx_base =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        for (unsigned int i = 0; i < chunk_size; ++i)
          {
            const size_type idx = idx_base + i * block_size;
            if (idx < N)
              val[idx] = s * val[idx] + a * V_val[idx] + b * W_val[idx];
          }
      }



      template <typename Number>
      __global__ void
      scale(Number *val, const Number *V_val, const size_type N)
      {
        const size_type idx_base =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        for (unsigned int i = 0; i < chunk_size; ++i)
          {
            const size_type idx = idx_base + i * block_size;
            if (idx < N)
              val[idx] *= V_val[idx];
          }
      }



      template <typename Number>
      __global__ void
      equ(Number *val, const Number a, const Number *V_val, const size_type N)
      {
        const size_type idx_base =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        for (unsigned int i = 0; i < chunk_size; ++i)
          {
            const size_type idx = idx_base + i * block_size;
            if (idx < N)
              val[idx] = a * V_val[idx];
          }
      }



      template <typename Number>
      __global__ void
      equ(Number *        val,
          const Number    a,
          const Number *  V_val,
          const Number    b,
          const Number *  W_val,
          const size_type N)
      {
        const size_type idx_base =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        for (unsigned int i = 0; i < chunk_size; ++i)
          {
            const size_type idx = idx_base + i * block_size;
            if (idx < N)
              val[idx] = a * V_val[idx] + b * W_val[idx];
          }
      }



      template <typename Number>
      __global__ void
      add_and_dot(Number *        res,
                  Number *        v1,
                  const Number *  v2,
                  const Number *  v3,
                  const Number    a,
                  const size_type N)
      {
        __shared__ Number res_buf[block_size];

        const unsigned int global_idx =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        const unsigned int local_idx = threadIdx.x;
        if (global_idx < N)
          {
            v1[global_idx] += a * v2[global_idx];
            res_buf[local_idx] =
              v1[global_idx] *
              Number(numbers::NumberTraits<Number>::conjugate(v3[global_idx]));
          }
        else
          res_buf[local_idx] = 0.;

        for (unsigned int i = 1; i < block_size; ++i)
          {
            const unsigned int idx = global_idx + i * block_size;
            if (idx < N)
              {
                v1[idx] += a * v2[idx];
                res_buf[local_idx] += v1[idx] * v3[idx];
              }
          }

        __syncthreads();

        reduce<Number, DotProduct<Number>>(
          res, res_buf, local_idx, global_idx, N);
      }



      template <typename Number>
      __global__ void
      set(Number *val, const Number s, const size_type N)
      {
        const size_type idx_base =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        for (unsigned int i = 0; i < chunk_size; ++i)
          {
            const size_type idx = idx_base + i * block_size;
            if (idx < N)
              val[idx] = s;
          }
      }



      template <typename Number>
      __global__ void
      set_permutated(Number *         val,
                     const Number *   v,
                     const size_type *indices,
                     const size_type  N)
      {
        const size_type idx_base =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        for (unsigned int i = 0; i < chunk_size; ++i)
          {
            const size_type idx = idx_base + i * block_size;
            if (idx < N)
              val[indices[idx]] = v[idx];
          }
      }



      template <typename Number>
      __global__ void
      add_permutated(Number *         val,
                     const Number *   v,
                     const size_type *indices,
                     const size_type  N)
      {
        const size_type idx_base =
          threadIdx.x + blockIdx.x * (blockDim.x * chunk_size);
        for (unsigned int i = 0; i < chunk_size; ++i)
          {
            const size_type idx = idx_base + i * block_size;
            if (idx < N)
              val[indices[idx]] += v[idx];
          }
      }
    } // namespace kernel
  }   // namespace CUDAWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
