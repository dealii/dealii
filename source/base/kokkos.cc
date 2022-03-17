// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2021 by the deal.II authors
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

namespace Impl{
 void cleanup_kokkos()
 {
    GrowingVectorMemory<
        LinearAlgebra::distributed::Vector<double, Kokkos::HostSpace>>::
        release_unused_memory();
      GrowingVectorMemory<
        LinearAlgebra::distributed::Vector<float, Kokkos::HostSpace>>::
        release_unused_memory();
#ifdef DEAL_II_WITH_CUDA
        GrowingVectorMemory<
        LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>>::
        release_unused_memory();
      GrowingVectorMemory<
        LinearAlgebra::distributed::Vector<float, MemorySpace::CUDA>>::
        release_unused_memory();
#endif
      Kokkos::finalize();
 }

 void ensure_kokkos_initialized()
 {
   if(!Kokkos::is_initialized())
	   GrowingVectorMemory<
        LinearAlgebra::distributed::Vector<double, Kokkos::HostSpace>>{};
      GrowingVectorMemory<
        LinearAlgebra::distributed::Vector<float, Kokkos::HostSpace>>{};
      #ifdef DEAL_II_WITH_CUDA
        GrowingVectorMemory<
        LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>>{};
      GrowingVectorMemory<
        LinearAlgebra::distributed::Vector<float, MemorySpace::CUDA>>{};
#endif
     Kokkos::initialize();
     std::atexit(cleanup_kokkos);
 }
}
DEAL_II_NAMESPACE_CLOSE


