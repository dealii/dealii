// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
   * that backend. Otherwise, the default memory space is the same as the
   * Host memory space.
   */
  struct Default
  {
    using kokkos_space = ::Kokkos::DefaultExecutionSpace::memory_space;
  };
} // namespace MemorySpace

DEAL_II_NAMESPACE_CLOSE

#endif
