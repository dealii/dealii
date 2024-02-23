// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_cuda_atomic_h
#define dealii_cuda_atomic_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CUDA

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace CUDAWrappers
  {
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
