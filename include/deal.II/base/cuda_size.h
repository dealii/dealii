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

#ifndef dealii_cuda_size_h
#define dealii_cuda_size_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  /**
   * Define the size of a block when launching a CUDA kernel. This number can be
   * changed depending on the architecture the code is running on.
   */
  constexpr int block_size = 512;

  /**
   * Define the size of chunk of data worked on by a thread. This number can be
   * changed depending on the architecture the code is running on.
   */
  constexpr int chunk_size = 1;

  /**
   * Define the number of threads in a warp.
   */
  constexpr int warp_size = 32;

  /**
   * Define the largest finite element degree that can be solved using
   * CUDAWrappers::MatrixFree. Changing this number will affect the amount of
   * constant memory being used.
   */
  constexpr unsigned int mf_max_elem_degree = 10;

  /**
   * Define the maximum number of valid CUDAWrappers::MatrixFree object.
   * Changing this number will affect the amount of constant memory being used.
   */
  constexpr unsigned int mf_n_concurrent_objects = 5;
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
