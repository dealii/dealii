// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#ifndef dealii_cuda_h
#define dealii_cuda_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/exceptions.h>

#ifdef DEAL_II_COMPILER_CUDA_AWARE
#  include <cusolverDn.h>
#  include <cusolverSp.h>
#  include <cusparse.h>

#  include <vector>

DEAL_II_NAMESPACE_OPEN
namespace Utilities
{
  /**
   * A namespace for utility structures for CUDA.
   */
  namespace CUDA
  {
    /**
     * Various CUDA APIs need an object to store internal data. This structure
     * creates, initializes, stores, and destroys these so-called handles for
     * the respective CUDA libraries used inside deal.II.
     */
    struct Handle
    {
      /**
       * Constructor. Initialize the handles for the different libraries.
       */
      Handle();

      /**
       * Copy constructor is deleted.
       */
      Handle(Handle const &) = delete;

      /**
       * Destructor. Destroy the handles.
       */
      ~Handle();

      /**
       * Pointer to an opaque cuSolverDN context.
       * The handle must be passed to every cuSolverDN library function.
       */
      cusolverDnHandle_t cusolver_dn_handle;

      /**
       * Pointer to an opaque cuSolverSP context.
       * The handle must be passed to every cuSolverSP library function.
       */
      cusolverSpHandle_t cusolver_sp_handle;

      /**
       * Pointer to an opaque cuSPARSE context.
       * The handle must be passed to every cuSPARSE library function.
       */
      cusparseHandle_t cusparse_handle;
    };

    /**
     * Allocate @p n_elements on the device.
     */
    template <typename T>
    inline void
    malloc(T *&pointer, const unsigned int n_elements)
    {
      cudaError_t cuda_error_code =
        cudaMalloc(&pointer, n_elements * sizeof(T));
      AssertCuda(cuda_error_code);
    }

    /**
     * Free memory on the device.
     */
    template <typename T>
    inline void
    free(T *&pointer)
    {
      cudaError_t cuda_error_code = cudaFree(pointer);
      AssertCuda(cuda_error_code);
      pointer = nullptr;
    }

    /**
     * Allocator to be used for `std::unique_ptr` pointing to device memory.
     */
    template <typename Number>
    Number *
    allocate_device_data(const std::size_t size)
    {
      Number *device_ptr;
      Utilities::CUDA::malloc(device_ptr, size);
      return device_ptr;
    }

    /**
     * Deleter to be used for `std::unique_ptr` pointing to device memory.
     */
    template <typename Number>
    void
    delete_device_data(Number *device_ptr) noexcept
    {
      const cudaError_t error_code = cudaFree(device_ptr);
      AssertNothrowCuda(error_code);
    }

    /**
     * Copy the device ArrayView @p in to the host ArrayView @p out.
     */
    template <typename T>
    inline void
    copy_to_host(const ArrayView<const T, MemorySpace::CUDA> &in,
                 ArrayView<T, MemorySpace::Host> &            out)
    {
      AssertDimension(in.size(), out.size());
      cudaError_t cuda_error_code = cudaMemcpy(out.data(),
                                               in.data(),
                                               in.size() * sizeof(T),
                                               cudaMemcpyDeviceToHost);
      AssertCuda(cuda_error_code);
    }

    /**
     * Copy the host ArrayView @p in to the device ArrayView @p out.
     */
    template <typename T>
    inline void
    copy_to_dev(const ArrayView<const T, MemorySpace::Host> &in,
                ArrayView<T, MemorySpace::CUDA> &            out)
    {
      AssertDimension(in.size(), out.size());
      cudaError_t cuda_error_code = cudaMemcpy(out.data(),
                                               in.data(),
                                               in.size() * sizeof(T),
                                               cudaMemcpyHostToDevice);
      AssertCuda(cuda_error_code);
    }

    /**
     * Copy the elements in @p pointer_dev to the host in @p vector_host.
     */
    template <typename T>
    inline void
    copy_to_host(const T *pointer_dev, std::vector<T> &vector_host)
    {
      ArrayView<const T, MemorySpace::CUDA> in(pointer_dev, vector_host.size());
      auto                                  out = make_array_view(vector_host);
      copy_to_host(in, out);
    }

    /**
     * Copy the elements in @p vector_host to the device in @p pointer_dev. The
     * memory needs to be allocate on the device before this function is called.
     */
    template <typename T>
    inline void
    copy_to_dev(const std::vector<T> &vector_host, T *pointer_dev)
    {
      auto                            in = make_array_view(vector_host);
      ArrayView<T, MemorySpace::CUDA> out(pointer_dev, vector_host.size());
      copy_to_dev(in, out);
    }
  } // namespace CUDA
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE
#endif
#endif
