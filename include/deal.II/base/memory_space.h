// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2023 by the deal.II authors
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


#ifndef dealii_memory_space_h
#define dealii_memory_space_h

#include <deal.II/base/config.h>

#include <Kokkos_Core.hpp>

DEAL_II_NAMESPACE_OPEN

/**
 */
namespace MemorySpace
{
  /**
   * Structure describing Host memory space.
   */
  struct Host
  {
    using kokkos_space = ::Kokkos::HostSpace;
  };

  /**
   * Structure describing the default memory space. If Kokkos was configured
   * with a GPU backend, the default memory space is the one corresponding to
   * that backend. Otherwise, the default memory space is the the same as the
   * Host memory space.
   */
  struct Default
  {
    using kokkos_space = ::Kokkos::DefaultExecutionSpace::memory_space;
  };

#ifdef DEAL_II_WITH_CUDA
  /**
   * Structure describing CUDA memory space.
   */
  using CUDA = Default;
#endif
} // namespace MemorySpace

DEAL_II_NAMESPACE_CLOSE

#endif
