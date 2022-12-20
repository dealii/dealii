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

#ifndef dealii_kokkos_h
#define dealii_kokkos_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * Records if Kokkos has been initialized by deal.II. The value stored is only
   * meaningful after ensure_kokkos_initialized() has been called.
   */
  extern bool dealii_initialized_kokkos;

  /**
   * Makes sure that Kokkos is initialized. Sets dealii_initialized_kokkos.
   */
  void
  ensure_kokkos_initialized();
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
