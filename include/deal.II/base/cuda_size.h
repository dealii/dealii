// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
