// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

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
#if DEAL_II_KOKKOS_VERSION_GTE(3, 7, 0)
        && !Kokkos::is_finalized()
#endif
    )
      {
        // only execute once
        static bool dummy = []() {
          dealii_initialized_kokkos = true;
#if DEAL_II_KOKKOS_VERSION_GTE(3, 7, 0)
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
