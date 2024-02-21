// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test that we can initialize and finalize Kokkos in user code.

#include <deal.II/base/kokkos.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_memory.h>

#include "../tests.h"

int
main()
{
  Kokkos::initialize();

  initlog();

  GrowingVectorMemory<
    LinearAlgebra::distributed::Vector<double, MemorySpace::Host>>{};
  GrowingVectorMemory<
    LinearAlgebra::distributed::Vector<float, MemorySpace::Host>>{};

  internal::ensure_kokkos_initialized();
  deallog << "Kokkos initialized by Kokkos: "
          << internal::dealii_initialized_kokkos << std::endl;

  Kokkos::finalize();

  return 0;
}
