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


#ifndef dealii_memory_space_h
#define dealii_memory_space_h

#include <deal.II/base/config.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <Kokkos_Core.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

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

  /**
   * Structure describing CUDA memory space.
   */
  // FIXME Only enable if CUDA is enabled in deal.II.
  using CUDA = Default;

} // namespace MemorySpace

DEAL_II_NAMESPACE_CLOSE

#endif
