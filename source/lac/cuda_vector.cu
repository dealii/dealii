// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/cuda_atomic.cuh>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/base/exceptions.h>
#include <cmath>

#ifdef DEAL_II_WITH_CUDA

DEAL_II_NAMESPACE_OPEN

#define BLOCK_SIZE 512
#define CHUNK_SIZE 8

namespace LinearAlgebra
{
  namespace CUDAWrappers
  {
    namespace internal
    {
      template <typename Number>
      __global__ void vec_scale(Number                                   *val,
                                const Number                              a,
                                const typename Vector<Number>::size_type  N)
      {
        const typename Vector<Number>::size_type idx_base = threadIdx.x +
                                                            blockIdx.x *
                                                            (blockDim.x*CHUNK_SIZE);
        for (unsigned int i=0; i<CHUNK_SIZE; ++i)
          {
            const typename Vector<Number>::size_type idx = idx_base +
                                                           i*BLOCK_SIZE;
            if (idx<N)
              val[idx] *= a;
          }
      }



      struct Binop_Addition
      {
        template <typename Number>
        __device__ static inline Number operation(const Number a,
                                                  const Number b)
        {
          return a+b;
        }
      };



      struct Binop_Subtraction
      {
        template <typename Number>
        __device__ static inline Number operation(const Number a,
                                                  const Number b)
        {
          return a-b;
        }
      };



      template <typename Number, typename Binop>
      __global__ void vector_bin_op(Number                                   *v1,
                                    Number                                   *v2,
                                    const typename Vector<Number>::size_type  N)
      {
        const typename Vector<Number>::size_type idx_base = threadIdx.x +
                                                            blockIdx.x *
                                                            (blockDim.x*CHUNK_SIZE);
        for (unsigned int i=0; i<CHUNK_SIZE; ++i)
          {
            const typename Vector<Number>::size_type idx = idx_base +
                                                           i*BLOCK_SIZE;
            if (idx<N)
              v1[idx] = Binop::operation(v1[idx],v2[idx]);
          }
      }



      template <typename Number>
      struct ElemSum
      {
        __device__ static Number reduction_op(const Number a, const Number b)
        {
          return (a + b);
        }

        __device__ static Number atomic_op(Number *dst, const Number a)
        {
          return atomicAdd_wrapper(dst, a);
        }

        __device__ static Number element_wise_op(const Number a)
        {
          return a;
        }

        __device__ static Number null_value()
        {
          return Number();
        }
      };



      template <typename Number>
      struct L1Norm
      {
        __device__ static Number reduction_op(const Number a, const Number b)
        {
          return (a + b);
        }

        __device__ static Number atomic_op(Number *dst, const Number a)
        {
          return atomicAdd_wrapper(dst, a);
        }

        __device__ static Number element_wise_op(const Number a)
        {
          return std::fabs(a);
        }

        __device__ static Number null_value()
        {
          return Number();
        }
      };



      template <typename Number>
      struct LInfty
      {
        __device__ static Number reduction_op(const Number a, const Number b)
        {
          if  (a > b)
            return a;
          else
            return b;
        }

        __device__ static Number atomic_op(Number *dst, const Number a)
        {
          return atomicMax_wrapper(dst, a);
        }

        __device__ static Number element_wise_op(const Number a)
        {
          return std::fabs(a);
        }

        __device__ static Number null_value()
        {
          return Number();
        }
      };



      template <typename Number, typename Operation>
      __device__ void reduce_within_warp(volatile Number                    *result_buffer,
                                         typename Vector<Number>::size_type  local_idx)
      {
        if (BLOCK_SIZE >= 64)
          result_buffer[local_idx] =
            Operation::reduction_op(result_buffer[local_idx],
                                    result_buffer[local_idx+32]);
        if (BLOCK_SIZE >= 32)
          result_buffer[local_idx] =
            Operation::reduction_op(result_buffer[local_idx],
                                    result_buffer[local_idx+16]);
        if (BLOCK_SIZE >= 16)
          result_buffer[local_idx] =
            Operation::reduction_op(result_buffer[local_idx],
                                    result_buffer[local_idx+8]);
        if (BLOCK_SIZE >= 8)
          result_buffer[local_idx] =
            Operation::reduction_op(result_buffer[local_idx],
                                    result_buffer[local_idx+4]);
        if (BLOCK_SIZE >= 4)
          result_buffer[local_idx] =
            Operation::reduction_op(result_buffer[local_idx],
                                    result_buffer[local_idx+2]);
        if (BLOCK_SIZE >= 2)
          result_buffer[local_idx] =
            Operation::reduction_op(result_buffer[local_idx],
                                    result_buffer[local_idx+1]);
      }



      template <typename Number, typename Operation>
      __device__ void reduce(Number                                   *result,
                             Number                                   *result_buffer,
                             const typename Vector<Number>::size_type  local_idx,
                             const typename Vector<Number>::size_type  global_idx,
                             const typename Vector<Number>::size_type  N)
      {
        for (typename Vector<Number>::size_type s=BLOCK_SIZE/2; s>32; s=s>>1)
          {
            if (local_idx < s)
              result_buffer[local_idx] = Operation::reduction_op(result_buffer[local_idx],
                                                                 result_buffer[local_idx+s]);
            __syncthreads();
          }

        if (local_idx < 32)
          reduce_within_warp<Number,Operation>(result_buffer, local_idx);

        if (local_idx == 0)
          Operation::atomic_op(result, result_buffer[0]);
      }



      template <typename Number, typename Operation>
      __global__ void reduction(Number       *result,
                                const Number *v,
                                const typename Vector<Number>::size_type N)
      {
        __shared__ Number result_buffer[BLOCK_SIZE];

        const typename Vector<Number>::size_type global_idx = threadIdx.x +
                                                              blockIdx.x*(blockDim.x*CHUNK_SIZE);
        const typename Vector<Number>::size_type local_idx = threadIdx.x;

        if (global_idx<N)
          result_buffer[local_idx] = Operation::element_wise_op(v[global_idx]);
        else
          result_buffer[local_idx] = Operation::null_value();

        __syncthreads();

        reduce<Number,Operation> (result, result_buffer, local_idx, global_idx, N);
      }



      template <typename Number>
      struct DotProduct
      {
        __device__ static Number binary_op(const Number a, const Number b)
        {
          return a*b;
        }

        __device__ static Number reduction_op(const Number a, const Number b)
        {
          return a+b;
        }

        __device__ static Number atomic_op(Number *dst, const Number a)
        {
          return atomicAdd_wrapper(dst, a);
        }

        __device__ static Number null_value()
        {
          return Number();
        }
      };



      template <typename Number, typename Operation>
      __global__ void double_vector_reduction(Number       *result,
                                              Number *v1,
                                              Number *v2,
                                              const typename Vector<Number>::size_type N)
      {
        __shared__ Number result_buffer[BLOCK_SIZE];

        const typename Vector<Number>::size_type global_idx = threadIdx.x +
                                                              blockIdx.x*(blockDim.x*CHUNK_SIZE);
        const typename Vector<Number>::size_type local_idx = threadIdx.x;

        if (global_idx<N)
          result_buffer[local_idx] = Operation::binary_op(v1[global_idx],v2[global_idx]);
        else
          result_buffer[local_idx] = Operation::null_value();

        for (unsigned int i=1; i<CHUNK_SIZE; ++i)
          {
            const typename Vector<Number>::size_type idx = global_idx +
                                                           i*BLOCK_SIZE;
            if (idx<N)
              result_buffer[local_idx] =
                Operation::reduction_op(result_buffer[local_idx],
                                        Operation::binary_op(v1[idx], v2[idx]));
          }

        __syncthreads();

        reduce<Number,Operation> (result,result_buffer,local_idx,global_idx,N);
      }



      template <typename Number>
      __global__ void vec_add(Number       *val,
                              const Number  a,
                              const typename Vector<Number>::size_type  N)
      {
        const typename Vector<Number>::size_type idx_base = threadIdx.x +
                                                            blockIdx.x *
                                                            (blockDim.x*CHUNK_SIZE);
        for (unsigned int i=0; i<CHUNK_SIZE; ++i)
          {
            const typename Vector<Number>::size_type idx = idx_base +
                                                           i*BLOCK_SIZE;
            if (idx<N)
              val[idx] += a;
          }
      }



      template <typename Number>
      __global__ void add_aV(Number       *val,
                             const Number  a,
                             Number       *V_val,
                             const typename Vector<Number>::size_type  N)
      {
        const typename Vector<Number>::size_type idx_base = threadIdx.x +
                                                            blockIdx.x *
                                                            (blockDim.x*CHUNK_SIZE);
        for (unsigned int i=0; i<CHUNK_SIZE; ++i)
          {
            const typename Vector<Number>::size_type idx = idx_base +
                                                           i*BLOCK_SIZE;
            if (idx<N)
              val[idx] += a*V_val[idx];
          }
      }



      template <typename Number>
      __global__ void add_aVbW(Number       *val,
                               const Number  a,
                               Number       *V_val,
                               const Number  b,
                               Number       *W_val,
                               const typename Vector<Number>::size_type  N)
      {
        const typename Vector<Number>::size_type idx_base = threadIdx.x +
                                                            blockIdx.x *
                                                            (blockDim.x*CHUNK_SIZE);
        for (unsigned int i=0; i<CHUNK_SIZE; ++i)
          {
            const typename Vector<Number>::size_type idx = idx_base +
                                                           i*BLOCK_SIZE;
            if (idx<N)
              val[idx] += a*V_val[idx] + b*W_val[idx];
          }
      }



      template <typename Number>
      __global__ void sadd(const Number  s,
                           Number       *val,
                           const Number  a,
                           const Number *V_val,
                           const typename Vector<Number>::size_type  N)
      {
        const typename Vector<Number>::size_type idx_base = threadIdx.x +
                                                            blockIdx.x *
                                                            (blockDim.x*CHUNK_SIZE);
        for (unsigned int i=0; i<CHUNK_SIZE; ++i)
          {
            const typename Vector<Number>::size_type idx = idx_base +
                                                           i*BLOCK_SIZE;
            if (idx<N)
              val[idx] = s*val[idx] + a*V_val[idx];
          }
      }



      template <typename Number>
      __global__ void scale(Number       *val,
                            const Number *V_val,
                            const typename Vector<Number>::size_type N)
      {
        const typename Vector<Number>::size_type idx_base = threadIdx.x +
                                                            blockIdx.x *
                                                            (blockDim.x*CHUNK_SIZE);
        for (unsigned int i=0; i<CHUNK_SIZE; ++i)
          {
            const typename Vector<Number>::size_type idx = idx_base +
                                                           i*BLOCK_SIZE;
            if (idx<N)
              val[idx] *= V_val[idx];
          }
      }



      template <typename Number>
      __global__ void equ(Number       *val,
                          const Number a,
                          const Number *V_val,
                          const typename Vector<Number>::size_type N)
      {
        const typename Vector<Number>::size_type idx_base = threadIdx.x +
                                                            blockIdx.x *
                                                            (blockDim.x*CHUNK_SIZE);
        for (unsigned int i=0; i<CHUNK_SIZE; ++i)
          {
            const typename Vector<Number>::size_type idx = idx_base +
                                                           i*BLOCK_SIZE;
            if (idx<N)
              val[idx] = a * V_val[idx];
          }
      }



      template <typename Number>
      __global__ void add_and_dot(Number       *res,
                                  Number       *v1,
                                  const Number *v2,
                                  const Number *v3,
                                  const Number  a,
                                  const typename Vector<Number>::size_type N)
      {
        __shared__ Number res_buf[BLOCK_SIZE];

        const unsigned int global_idx = threadIdx.x + blockIdx.x *
                                        (blockDim.x*CHUNK_SIZE);
        const unsigned int local_idx = threadIdx.x;
        if (global_idx < N)
          {
            v1[global_idx] += a*v2[global_idx];
            res_buf[local_idx] = v1[global_idx]*v3[global_idx];
          }
        else
          res_buf[local_idx] = 0.;

        for (unsigned int i=1; i<BLOCK_SIZE; ++i)
          {
            const unsigned int idx = global_idx + i*BLOCK_SIZE;
            if (idx < N)
              {
                v1[idx] += a*v2[idx];
                res_buf[local_idx] += v1[idx]*v3[idx];
              }
          }

        __syncthreads();

        reduce<Number, DotProduct<Number>> (res, res_buf, local_idx,
                                            global_idx, N);
      }
    }



    template <typename Number>
    Vector<Number>::Vector()
      :
      val(nullptr),
      n_elements(0)
    {}



    template <typename Number>
    Vector<Number>::Vector(const Vector<Number> &V)
      :
      n_elements(V.n_elements)
    {
      // Allocate the memory
      cudaError_t error_code = cudaMalloc(&val, n_elements*sizeof(Number));
      AssertCuda(error_code);
      // Copy the values.
      error_code = cudaMemcpy(val, V.val,n_elements*sizeof(Number),
                              cudaMemcpyDeviceToDevice);
      AssertCuda(error_code);
    }



    template <typename Number>
    Vector<Number>::Vector(const size_type n)
      :
      n_elements(n)
    {
      // Allocate the memory
      cudaError_t error_code = cudaMalloc(&val, n_elements*sizeof(Number));
      AssertCuda(error_code);
    }



    template <typename Number>
    Vector<Number>::~Vector()
    {
      if (val != nullptr)
        {
          cudaError_t error_code = cudaFree(val);
          AssertCuda(error_code);
          val = nullptr;
          n_elements = 0;
        }
    }



    template <typename Number>
    void Vector<Number>::reinit(const size_type n,
                                const bool      omit_zeroing_entries)
    {
      // Resize the underlying array if necessary
      if (n == 0)
        {
          if (val != nullptr)
            {
              cudaError_t error_code = cudaFree(val);
              AssertCuda(error_code);
              val = nullptr;
            }
        }
      else
        {
          if (n_elements != n)
            {
              cudaError_t error_code = cudaFree(val);
              AssertCuda(error_code);
            }

          cudaError_t error_code = cudaMalloc(&val, n*sizeof(Number));
          AssertCuda(error_code);

          // If necessary set the elements to zero
          if (omit_zeroing_entries == false)
            {
              cudaError_t error_code = cudaMemset(val, 0,
                                                  n_elements*sizeof(Number));
              AssertCuda(error_code);
            }
        }
      n_elements = n;
    }



    template <typename Number>
    void Vector<Number>::reinit(const VectorSpaceVector<Number> &V,
                                const bool omit_zeroing_entries)
    {
      reinit(V.size(), omit_zeroing_entries);
    }



    template <typename Number>
    void Vector<Number>::import(const ReadWriteVector<Number> &V,
                                VectorOperation::values operation,
                                std::shared_ptr<const CommunicationPatternBase> )
    {
      if (operation == VectorOperation::insert)
        {
          cudaError_t error_code = cudaMemcpy(val, V.begin(),
                                              n_elements*sizeof(Number),
                                              cudaMemcpyHostToDevice);
          AssertCuda(error_code);
        }
      else
        {
          // Create a temporary vector on the device
          Number *tmp;
          cudaError_t error_code = cudaMalloc(&tmp, n_elements*sizeof(Number));
          AssertCuda(error_code);

          // Copy the vector from the host to the temporary vector on the device
          error_code = cudaMemcpy(&tmp[0], V.begin(), n_elements*sizeof(Number),
                                  cudaMemcpyHostToDevice);
          AssertCuda(error_code);

          // Add the two vectors
          const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);

          internal::vector_bin_op<Number,internal::Binop_Addition>
          <<<n_blocks,BLOCK_SIZE>>>(val, tmp, n_elements);
          // Check that the kernel was launched correctly
          AssertCuda(cudaGetLastError());
          // Check that there was no problem during the execution of the kernel
          AssertCuda(cudaDeviceSynchronize());

          // Delete the temporary vector
          error_code = cudaFree(tmp);
          AssertCuda(error_code);
        }
    }



    template <typename Number>
    Vector<Number> &Vector<Number>::operator= (const Number s)
    {
      Assert(s == Number(), ExcMessage("Onlyt 0 can be assigned to a vector."));
      (void)s;

      cudaError_t error_code = cudaMemset(val, 0, n_elements*sizeof(Number));
      AssertCuda(error_code);

      return *this;
    }



    template <typename Number>
    Vector<Number> &Vector<Number>::operator*= (const Number factor)
    {
      AssertIsFinite(factor);
      const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);
      internal::vec_scale<Number> <<<n_blocks,BLOCK_SIZE>>>(val,
                                                            factor, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());

      return *this;
    }



    template <typename Number>
    Vector<Number> &Vector<Number>::operator/= (const Number factor)
    {
      AssertIsFinite(factor);
      Assert(factor!=Number(0.), ExcZero());
      const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);
      internal::vec_scale<Number> <<<n_blocks,BLOCK_SIZE>>>(val,
                                                            1./factor, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());

      return *this;
    }



    template <typename Number>
    Vector<Number> &Vector<Number>::operator+= (const VectorSpaceVector<Number> &V)
    {
      // Check that casting will work
      Assert(dynamic_cast<const Vector<Number>*>(&V)!=nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If it fails, it throw an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
      Assert(down_V.size()==this->size(),
             ExcMessage("Cannot add two vectors with different numbers of elements"));

      const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);

      internal::vector_bin_op<Number,internal::Binop_Addition>
      <<<n_blocks,BLOCK_SIZE>>>(val, down_V.val, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());

      return *this;
    }



    template <typename Number>
    Vector<Number> &Vector<Number>::operator-= (const VectorSpaceVector<Number> &V)
    {
      // Check that casting will work
      Assert(dynamic_cast<const Vector<Number>*>(&V)!=nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
      Assert(down_V.size()==this->size(),
             ExcMessage("Cannot add two vectors with different numbers of elements."));

      const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);

      internal::vector_bin_op<Number,internal::Binop_Subtraction>
      <<<n_blocks,BLOCK_SIZE>>>(val, down_V.val, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());

      return *this;
    }



    template <typename Number>
    Number Vector<Number>::operator* (const VectorSpaceVector<Number> &V) const
    {
      // Check that casting will work
      Assert(dynamic_cast<const Vector<Number>*>(&V)!=nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
      Assert(down_V.size()==this->size(),
             ExcMessage("Cannot add two vectors with different numbers of elements"));

      Number *result_device;
      cudaError_t error_code = cudaMalloc(&result_device, n_elements*sizeof(Number));
      AssertCuda(error_code);
      error_code = cudaMemset(result_device, Number(), sizeof(Number));

      const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);
      internal::double_vector_reduction<Number, internal::DotProduct<Number>>
          <<<dim3(n_blocks,1),dim3(BLOCK_SIZE)>>> (result_device, val,
                                                   down_V.val,
                                                   static_cast<unsigned int>(n_elements));

      // Copy the result back to the host
      Number result;
      error_code = cudaMemcpy(&result, result_device, sizeof(Number),
                              cudaMemcpyDeviceToHost);
      AssertCuda(error_code);
      // Free the memory on the device
      error_code = cudaFree(result_device);
      AssertCuda(error_code);

      return result;
    }



    template <typename Number>
    void Vector<Number>::add(const Number a)
    {
      AssertIsFinite(a);
      const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);
      internal::vec_add<Number> <<<n_blocks,BLOCK_SIZE>>>(val, a,
                                                          n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());
    }



    template <typename Number>
    void Vector<Number>::add(const Number a, const VectorSpaceVector<Number> &V)
    {
      AssertIsFinite(a);

      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number>*>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throw an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
      Assert(down_V.size() == this->size(),
             ExcMessage("Cannot add two vectors with different numbers of elements."));

      const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);
      internal::add_aV<Number> <<<dim3(n_blocks,1),dim3(BLOCK_SIZE)>>> (val,
          a, down_V.val, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());
    }



    template <typename Number>
    void Vector<Number>::add(const Number a, const VectorSpaceVector<Number> &V,
                             const Number b, const VectorSpaceVector<Number> &W)
    {
      AssertIsFinite(a);
      AssertIsFinite(b);

      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number>*>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throw an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
      Assert(down_V.size() == this->size(),
             ExcMessage("Cannot add two vectors with different numbers of elements."));

      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number>*>(&W) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throw an exception.
      const Vector<Number> &down_W = dynamic_cast<const Vector<Number>&>(W);
      Assert(down_W.size() == this->size(),
             ExcMessage("Cannot add two vectors with different numbers of elements."));

      const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);
      internal::add_aVbW<Number> <<<dim3(n_blocks,1),dim3(BLOCK_SIZE)>>> (val,
          a, down_V.val, b, down_W.val, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());
    }



    template <typename Number>
    void Vector<Number>::sadd(const Number s, const Number a,
                              const VectorSpaceVector<Number> &V)
    {
      AssertIsFinite(s);
      AssertIsFinite(a);

      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number>*>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throw an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
      Assert(down_V.size() == this->size(),
             ExcMessage("Cannot add two vectors with different numbers of elements."));

      const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);
      internal::sadd<Number> <<<dim3(n_blocks,1),dim3(BLOCK_SIZE)>>> (s, val,
          a, down_V.val, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());
    }



    template <typename Number>
    void Vector<Number>::scale(const VectorSpaceVector<Number> &scaling_factors)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number>*>(&scaling_factors) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throw an exception.
      const Vector<Number> &down_scaling_factors =
        dynamic_cast<const Vector<Number>&>(scaling_factors);
      Assert(down_scaling_factors.size() == this->size(),
             ExcMessage("Cannot scale two vectors with different numbers of elements."));

      const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);
      internal::scale<Number> <<<dim3(n_blocks,1),dim3(BLOCK_SIZE)>>> (val,
          down_scaling_factors.val, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());
    }



    template <typename Number>
    void Vector<Number>::equ(const Number a, const VectorSpaceVector<Number> &V)
    {
      AssertIsFinite(a);

      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number>*>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throw an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
      Assert(down_V.size() == this->size(),
             ExcMessage("Cannot assign two vectors with different numbers of elements."));

      const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);
      internal::equ<Number> <<<dim3(n_blocks,1),dim3(BLOCK_SIZE)>>> (val, a,
          down_V.val, n_elements);

      // Check that the kernel was launched correctly
      AssertCuda(cudaGetLastError());
      // Check that there was no problem during the execution of the kernel
      AssertCuda(cudaDeviceSynchronize());
    }



    template <typename Number>
    typename Vector<Number>::value_type Vector<Number>::mean_value() const
    {
      Number *result_device;
      cudaError_t error_code = cudaMalloc(&result_device, sizeof(Number));
      AssertCuda(error_code);
      error_code = cudaMemset(result_device, Number(), sizeof(Number));

      const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);
      internal::reduction<Number, internal::ElemSum<Number>>
                                                          <<<dim3(n_blocks,1),dim3(BLOCK_SIZE)>>> (
                                                            result_device, val,
                                                            n_elements);

      // Copy the result back to the host
      Number result;
      error_code = cudaMemcpy(&result, result_device, sizeof(Number),
                              cudaMemcpyDeviceToHost);
      AssertCuda(error_code);
      // Free the memory on the device
      error_code = cudaFree(result_device);
      AssertCuda(error_code);

      return result/static_cast<typename Vector<Number>::value_type>(n_elements);
    }



    template <typename Number>
    typename Vector<Number>::real_type Vector<Number>::l1_norm() const
    {
      Number *result_device;
      cudaError_t error_code = cudaMalloc(&result_device, sizeof(Number));
      AssertCuda(error_code);
      error_code = cudaMemset(result_device, Number(), sizeof(Number));

      const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);
      internal::reduction<Number, internal::L1Norm<Number>>
                                                         <<<dim3(n_blocks,1),dim3(BLOCK_SIZE)>>> (
                                                           result_device, val,
                                                           n_elements);

      // Copy the result back to the host
      Number result;
      error_code = cudaMemcpy(&result, result_device, sizeof(Number),
                              cudaMemcpyDeviceToHost);
      AssertCuda(error_code);
      // Free the memory on the device
      error_code = cudaFree(result_device);
      AssertCuda(error_code);

      return result;
    }



    template <typename Number>
    typename Vector<Number>::real_type Vector<Number>::l2_norm() const
    {
      return std::sqrt((*this)*(*this));
    }



    template <typename Number>
    typename Vector<Number>::real_type Vector<Number>::linfty_norm() const
    {
      Number *result_device;
      cudaError_t error_code = cudaMalloc(&result_device, sizeof(Number));
      AssertCuda(error_code);
      error_code = cudaMemset(result_device, Number(), sizeof(Number));

      const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);
      internal::reduction<Number, internal::LInfty<Number>>
                                                         <<<dim3(n_blocks,1),dim3(BLOCK_SIZE)>>> (
                                                           result_device, val,
                                                           n_elements);

      // Copy the result back to the host
      Number result;
      error_code = cudaMemcpy(&result, result_device, sizeof(Number),
                              cudaMemcpyDeviceToHost);
      AssertCuda(error_code);
      // Free the memory on the device
      error_code = cudaFree(result_device);
      AssertCuda(error_code);

      return result;
    }



    template <typename Number>
    Number Vector<Number>::add_and_dot(const Number                     a,
                                       const VectorSpaceVector<Number> &V,
                                       const VectorSpaceVector<Number> &W)
    {
      AssertIsFinite(a);

      // Check that casting will work
      Assert(dynamic_cast<const Vector<Number>*>(&V)!=nullptr,
             ExcVectorTypeNotCompatible());
      Assert(dynamic_cast<const Vector<Number>*>(&W)!=nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V and W. If it fails, throw an exceptiion.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
      Assert(down_V.size() == this->size(),
             ExcMessage("Vector V has the wrong size."));
      const Vector<Number> &down_W = dynamic_cast<const Vector<Number>&>(W);
      Assert(down_W.size() == this->size(),
             ExcMessage("Vector W has the wrong size."));

      Number *res_d;
      cudaError_t error_code = cudaMalloc(&res_d, sizeof(Number));
      AssertCuda(error_code);
      error_code = cudaMemset(res_d, 0., sizeof(Number));
      AssertCuda(error_code);

      const int n_blocks = 1 + (n_elements-1)/(CHUNK_SIZE*BLOCK_SIZE);
      internal::add_and_dot<Number> <<<dim3(n_blocks,1),dim3(BLOCK_SIZE)>>>(
        res_d, val, down_V.val, down_W.val, a, n_elements);

      Number res;
      error_code = cudaMemcpy(&res, res_d, sizeof(Number), cudaMemcpyDeviceToHost);
      AssertCuda(error_code);
      error_code = cudaFree(res_d);

      return res;
    }



    template <typename Number>
    void Vector<Number>::print(std::ostream       &out,
                               const unsigned int  precision,
                               const bool          scientific,
                               const bool          ) const
    {
      AssertThrow(out, ExcIO());
      std::ios::fmtflags old_flags = out.flags();
      unsigned int old_precision = out.precision (precision);

      out.precision (precision);
      if (scientific)
        out.setf (std::ios::scientific, std::ios::floatfield);
      else
        out.setf (std::ios::fixed, std::ios::floatfield);

      out << "IndexSet: ";
      complete_index_set(n_elements).print(out);
      out << std::endl;

      // Copy the vector to the host
      Number *cpu_val = new Number[n_elements];
      cudaError_t error_code = cudaMemcpy(cpu_val, val,
                                          n_elements*sizeof(Number),
                                          cudaMemcpyHostToDevice);
      AssertCuda(error_code);
      for (unsigned int i=0; i<n_elements; ++i)
        out << cpu_val[i] << std::endl;
      out << std::flush;
      delete [] cpu_val;
      cpu_val = nullptr;


      AssertThrow (out, ExcIO());
      // reset output format
      out.flags (old_flags);
      out.precision(old_precision);
    }



    template <typename Number>
    std::size_t Vector<Number>::memory_consumption() const
    {
      std::size_t memory = sizeof(*this);
      memory += sizeof (Number) * static_cast<std::size_t>(n_elements);

      return memory;
    }



    // Explicit Instanationation
    template class Vector<float>;
    template class Vector<double>;
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
