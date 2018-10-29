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

#ifndef dealii_cuda_size_h
#define dealii_cuda_size_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  /**
   * Define the size of a block when launching a CUDA kernel.
   */
  constexpr int block_size = 512;

  /**
   * Define the size of chunk of data worked on by a thread.
   */
  constexpr int chunk_size = 8;
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
