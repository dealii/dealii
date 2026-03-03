// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2016 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_distributed_amr_h
#define dealii_distributed_amr_h

#ifdef DEAL_II_WITH_P4EST
#  include <deal.II/distributed/p4est_wrappers.h>
#endif // DEAL_II_WITH_P4EST

namespace dealii
{
  namespace internal
  {
    namespace amr
    {
#if defined(DEAL_II_WITH_P4EST)
      using namespace dealii::internal::p4est;
#else
#  error DEAL_II_WITH_P4EST required
#endif
    } // namespace amr
  }   // namespace internal
} // namespace dealii


#endif // dealii_distributed_amr_h
