// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

#ifndef dealii_cuda_atomic_h
#define dealii_cuda_atomic_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_COMPILER_CUDA_AWARE

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace CUDAWrappers
  {
    /**
     * Provide atomicAdd for floats.
     *
     * @deprecated Use atomicAdd(address, val) directly.
     * @ingroup CUDAWrappers
     */
    DEAL_II_DEPRECATED
    inline __device__ float
    atomicAdd_wrapper(float *address, float val)
    {
      return atomicAdd(address, val);
    }



    /**
     * Provide atomicAdd for doubles.
     *
     * @deprecated Use atomicAdd(address, val) directly.
     * @ingroup CUDAWrappers
     */
    DEAL_II_DEPRECATED
    inline __device__ double
    atomicAdd_wrapper(double *address, double val)
    {
      return atomicAdd(address, val);
    }



    /**
     * Provide atomicMax for floats.
     *
     * @ingroup CUDAWrappers
     */
    inline __device__ float
    atomicMax_wrapper(float *address, float val)
    {
      int *address_as_int = reinterpret_cast<int *>(address);
      int  old            = *address_as_int, assumed;
      do
        {
          assumed = old;
          old     = atomicCAS(address_as_int,
                          assumed,
                          atomicMax(address_as_int, __float_as_int(val)));
        }
      while (assumed != old);

      return __longlong_as_double(old);
    }



    /**
     * Provide atomicMax for doubles.
     *
     * @ingroup CUDAWrappers
     */
    inline __device__ double
    atomicMax_wrapper(double *address, double val)
    {
      unsigned long long int *address_as_ull =
        reinterpret_cast<unsigned long long int *>(address);
      unsigned long long int old = *address_as_ull, assumed;
      do
        {
          assumed = old;
          old     = atomicCAS(address_as_ull,
                          assumed,
                          atomicMax(address_as_ull,
                                    static_cast<unsigned long long int>(
                                      __double_as_longlong(val))));
        }
      while (assumed != old);

      return __longlong_as_double(old);
    }
  } // namespace CUDAWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
