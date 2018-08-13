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

#ifndef dealii_cuda_h
#define dealii_cuda_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#if defined(DEAL_II_WITH_CUDA) && defined(__CUDACC__)
#  include <cusolverDn.h>
#  include <cusolverSp.h>
#  include <cusparse.h>
#endif

#include <vector>

DEAL_II_NAMESPACE_OPEN
namespace Utilities
{
  /**
   * A namespace for utility structures for CUDA.
   */
  namespace CUDA
  {
    /**
     * This structure creates, stores, and destroys the handles of the different
     * CUDA libraries used inside deal.II.
     */
    struct Handle
    {
      /**
       * Constructor. Create the handles for the different libraries.
       */
      Handle();

      /**
       * Copy constructor is deleted.
       */
      Handle(Handle const &) = delete;

      /**
       * Destructor. Destroy the handles and free all the memory allocated by
       * GrowingVectorMemory.
       */
      ~Handle();

#if defined(DEAL_II_WITH_CUDA) && defined(__CUDACC__)
      cusolverDnHandle_t cusolver_dn_handle;

      cusolverSpHandle_t cusolver_sp_handle;

      /**
       * Handle to the cuSPARSE library.
       */
      cusparseHandle_t cusparse_handle;
#endif
    };

    /**
     * Allocate @p n_elements on the device.
     */
    template <typename T>
    inline void
    malloc(T *&pointer, const unsigned int n_elements)
    {
#if defined(DEAL_II_WITH_CUDA) && defined(__CUDACC__)
      cudaError_t cuda_error_code =
        cudaMalloc(&pointer, n_elements * sizeof(T));
      AssertCuda(cuda_error_code);
#else
      (void)pointer;
      (void)n_elements;
#endif
    }

    /**
     * Free memory on the device.
     */
    template <typename T>
    inline void
    free(T *&pointer)
    {
#if defined(DEAL_II_WITH_CUDA) && defined(__CUDACC__)
      cudaError_t cuda_error_code = cudaFree(pointer);
      AssertCuda(cuda_error_code);
      pointer = nullptr;
#else
      (void)pointer;
#endif
    }

    /**
     * Copy the elements in @p pointer_dev to the host in @p vector_host.
     */
    template <typename T>
    inline void
    copy_to_host(const T *pointer_dev, std::vector<T> &vector_host)
    {
#if defined(DEAL_II_WITH_CUDA) && defined(__CUDACC__)
      cudaError_t cuda_error_code = cudaMemcpy(vector_host.data(),
                                               pointer_dev,
                                               vector_host.size() * sizeof(T),
                                               cudaMemcpyDeviceToHost);
      AssertCuda(cuda_error_code);
#else
      (void)pointer_dev;
      (void)vector_host;
#endif
    }

    /**
     * Copy the elements in @p vector_host to the device in @p pointer_dev. The
     * memory needs to be allocate on the device before this function is called.
     */
    template <typename T>
    inline void
    copy_to_dev(const std::vector<T> &vector_host, T *pointer_dev)
    {
#if defined(DEAL_II_WITH_CUDA) && defined(__CUDACC__)
      cudaError_t cuda_error_code = cudaMemcpy(pointer_dev,
                                               vector_host.data(),
                                               vector_host.size() * sizeof(T),
                                               cudaMemcpyHostToDevice);
      AssertCuda(cuda_error_code);
#else
      (void)vector_host;
      (void)pointer_dev;
#endif
    }
  } // namespace CUDA
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif
