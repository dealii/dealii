// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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

#ifndef dealii_memory_space_utils_h
#define dealii_memory_space_utils_h

#include <deal.II/base/config.h>

#include <deal.II/base/cuda.h>
#include <deal.II/base/cuda_size.h>
#include <deal.II/base/memory_space.h>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  /**
   * A namespace for utility functions that are independent of the memory space.
   *
   * @ingroup utilities
   */
  namespace MemorySpace
  {
#ifdef DEAL_II_COMPILER_CUDA_AWARE
    /**
     * Allocate @p size elements of type `Number` on the device.
     */
    template <typename Number>
    Number *
    allocate_data(const std::size_t size, const ::dealii::MemorySpace::CUDA &)
    {
      return ::dealii::Utilities::CUDA::allocate_device_data<Number>(size);
    }

    /**
     * Free the pointer.
     */
    template <typename Number>
    void
    delete_data(Number *ptr, const ::dealii::MemorySpace::CUDA &) noexcept
    {
      ::dealii::Utilities::CUDA::delete_device_data(ptr);
    }

    /**
     * Copy @p size elements fom @p in on the device to @p out on the host.
     */
    template <typename Number>
    void
    copy(const Number *in,
         const ::dealii::MemorySpace::CUDA &,
         Number *out,
         const ::dealii::MemorySpace::Host &,
         std::size_t size)
    {
      ArrayView<const Number, ::dealii::MemorySpace::CUDA> in_view(in, size);
      ArrayView<Number, ::dealii::MemorySpace::Host>       out_view(out, size);
      ::dealii::Utilities::CUDA::copy_to_host(in_view, out_view);
    }

    /**
     * Copy @p size elements fom @p in on the host to @p out on the device.
     */
    template <typename Number>
    void
    copy(const Number *in,
         const ::dealii::MemorySpace::Host &,
         Number *out,
         const ::dealii::MemorySpace::CUDA &,
         const std::size_t size)
    {
      ArrayView<const Number, ::dealii::MemorySpace::Host> in_view(in, size);
      ArrayView<Number, ::dealii::MemorySpace::CUDA>       out_view(out, size);
      ::dealii::Utilities::CUDA::copy_to_dev(in_view, out_view);
    }

    /**
     * Copy @p size elements fom @p in on the device to @p out on the device.
     */
    template <typename Number>
    void
    copy(const Number *in,
         const ::dealii::MemorySpace::CUDA &,
         Number *out,
         const ::dealii::MemorySpace::CUDA &,
         const std::size_t size)
    {
      cudaError_t cuda_error_code =
        cudaMemcpy(out, in, size * sizeof(Number), cudaMemcpyDeviceToDevice);
      AssertCuda(cuda_error_code);
    }

    namespace internal
    {
      /**
       * Helper function.
       */
      template <typename Functor>
      __global__ void
      for_each_impl(const unsigned int size, Functor f)
      {
        unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
        if (i < size)
          {
            f(i);
          }
      }
    } // namespace internal

    /**
     * Apply the functor @p f to the range [0,size). The first argument of the
     * function is an executor. Here the executor is the CUDA memory space and
     * the loop is executed by launching a kernel. In case the Functor @p f is
     * a lambda function, the code should look as follows
     * ```
     * for_each_index(MemorySpace::CUDA{}, size, [=] DEAL_II_CUDA_HOST_DEV (int
     * i)
     * {....});
     * ```
     * There is a few things to note:
     *  1. The capture has to be done by value, capture by reference will not
     *  work
     *  2. DEAL_II_CUDA_HOST_DEV expands to __host__ __device__ when nvcc is the
     *  compiler, otherwise it expends to nothing. Therefore, it is safe to add
     *  the macro even if the memory space is Host{}
     *  3. To enable lambda support in CUDA, you need to compile deal.II with
     *  --expt-extended-lambda or --extended-lambda depending of the version of
     *  CUDA that you are using
     *  4. Lambda function in CUDA have several limitations including:
     *    - they cannot be called in a constructor or in a private function
     *    - they cannot use private data member
     */
    template <typename Functor>
    void
    for_each_index(const dealii::MemorySpace::CUDA &,
                   const unsigned int size,
                   Functor            f)
    {
      const int n_blocks = 1 + size / dealii::CUDAWrappers::block_size;
      internal::for_each_impl<<<n_blocks, dealii::CUDAWrappers::block_size>>>(
        size, f);
      AssertCudaKernel();
    }
#endif

    /**
     * Allocate @p size elements of type `Number` on the host.
     */
    template <typename Number>
    Number *
    allocate_data(const std::size_t size, const ::dealii::MemorySpace::Host &)
    {
      Number *host_ptr = new Number[size];

      return host_ptr;
    }

    /**
     * Free the pointer.
     */
    template <typename Number>
    void
    delete_data(Number *ptr, const ::dealii::MemorySpace::Host &) noexcept
    {
      delete[] ptr;
    }

    /**
     * Copy @p size elements fom @p in on the host to @p out on the host.
     */
    template <typename Number>
    void
    copy(const Number *in,
         const ::dealii::MemorySpace::Host &,
         Number *out,
         const ::dealii::MemorySpace::Host &,
         const std::size_t size)
    {
      std::memcpy(out, in, size * sizeof(Number));
    }

    /**
     * Apply the functor @p f to the range [0,size). The first argument of the
     * function is an executor. Here the executor is the Host memory space and
     * the loop is executed in serial. This function accepts a lambda function
     * or a functor.
     */
    template <typename Functor>
    inline void
    for_each_index(const dealii::MemorySpace::Host &,
                   const unsigned int size,
                   Functor            f)
    {
      for (unsigned int i = 0; i < size; ++i)
        f(i);
    }

    /**
     * Allocate memory.
     */
    template <typename Number, typename MemorySpaceType>
    Number *
    allocate_data(const std::size_t size)
    {
      return allocate_data<Number>(size, MemorySpaceType{});
    }

    /**
     * Free the pointer.
     */
    template <typename Number, typename MemorySpaceType>
    void
    delete_data(Number *ptr) noexcept
    {
      delete_data<Number>(ptr, MemorySpaceType{});
    }

    /**
     * Copy @p in to @p out.
     */
    template <typename Number,
              typename MemorySpaceType,
              typename MemorySpaceType2>
    void
    copy(const ArrayView<const Number, MemorySpaceType> &in,
         ArrayView<Number, MemorySpaceType2> &           out)
    {
      AssertDimension(in.size(), out.size());
      copy(in.data(),
           MemorySpaceType{},
           out.data(),
           MemorySpaceType2{},
           in.size());
    }

  } // namespace MemorySpace
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif
