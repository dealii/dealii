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
