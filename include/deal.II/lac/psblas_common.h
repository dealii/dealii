// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_psblas_common_h
#define dealii_psblas_common_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PSBLAS

#  include <psb_base_cbind.h>
#  include <psb_c_base.h>
#  include <psb_c_dbase.h>

DEAL_II_NAMESPACE_OPEN

namespace PSCToolkitWrappers
{

  namespace internal
  {
    /*
     * Custom deleter for PSBLAS descriptor.
     */
    struct PSBLASDescriptorDeleter
    {
      void
      operator()(psb_c_descriptor *p) const
      {
        if (p)
          psb_c_cdfree(p);
      }
    };

    /**
     * Enum to indicate the state of the vector (building or assembled).
     */

    enum State
    {
      /**
       * State entered after the default constructor, before any allocation.
       * In this state, no operations are possible.
       */
      Default,
      /**
       * State entered after the first allocation, and before the first
       * assembly; in this state it is possible to add communication
       * requirements among different processes.
       */
      Build,
      /*
       * State entered after the assembly; computations such as matrix-vector
       * products, are only possible in this state.
       */
      Assembled
    };

  } // namespace internal
} // namespace PSCToolkitWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PSBLAS
#endif