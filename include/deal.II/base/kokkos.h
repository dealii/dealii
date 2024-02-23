// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
