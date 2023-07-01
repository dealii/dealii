// ---------------------------------------------------------------------
//
// Copyright (C) 2022 - 2023 by the deal.II authors
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
#include <deal.II/base/multithread_info.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_memory.h>

#include <Kokkos_Core.hpp>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  bool dealii_initialized_kokkos = false;

  void
  ensure_kokkos_initialized()
  {
    if (!Kokkos::is_initialized()
#if KOKKOS_VERSION >= 30700
        && !Kokkos::is_finalized()
#endif
    )
      {
        // only execute once
        static bool dummy = [] {
          dealii_initialized_kokkos = true;
#if KOKKOS_VERSION >= 30700
          const auto settings =
            Kokkos::InitializationSettings().set_num_threads(
              MultithreadInfo::n_threads());
#else
          const Kokkos::InitArguments settings(MultithreadInfo::n_threads());
#endif
          Kokkos::initialize(settings);
          std::atexit(Kokkos::finalize);
          return true;
        }();
        (void)dummy;
      }
  }
} // namespace internal
DEAL_II_NAMESPACE_CLOSE
