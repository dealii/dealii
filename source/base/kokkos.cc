// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

#include <deal.II/base/kokkos.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_memory.h>

#include <Kokkos_Core.hpp>

DEAL_II_NAMESPACE_OPEN

namespace Impl
{
  void
  ensure_kokkos_initialized()
  {
    if (!Kokkos::is_initialized())
      {
        GrowingVectorMemory<
          LinearAlgebra::distributed::Vector<double, MemorySpace::Host>>{};
        GrowingVectorMemory<
          LinearAlgebra::distributed::Vector<float, MemorySpace::Host>>{};
        GrowingVectorMemory<
          LinearAlgebra::distributed::Vector<double, MemorySpace::Device>>{};
        GrowingVectorMemory<
          LinearAlgebra::distributed::Vector<float, MemorySpace::Device>>{};
        Kokkos::push_finalize_hook(
          GrowingVectorMemory<
            LinearAlgebra::distributed::Vector<double, MemorySpace::Host>>::
            release_unused_memory);
        Kokkos::push_finalize_hook(
          GrowingVectorMemory<
            LinearAlgebra::distributed::Vector<float, MemorySpace::Host>>::
            release_unused_memory);
        Kokkos::push_finalize_hook(
          GrowingVectorMemory<
            LinearAlgebra::distributed::Vector<double, MemorySpace::Device>>::
            release_unused_memory);
        Kokkos::push_finalize_hook(
          GrowingVectorMemory<
            LinearAlgebra::distributed::Vector<float, MemorySpace::Device>>::
            release_unused_memory);
        Kokkos::initialize();
        std::atexit(Kokkos::finalize);
      }
  }
} // namespace Impl
DEAL_II_NAMESPACE_CLOSE
