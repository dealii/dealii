// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef amr_h
#define amr_h

#ifdef DEAL_II_WITH_P4EST
#  include <deal.II/distributed/p4est_wrappers.h>
#endif // DEAL_II_WITH_P4EST
#ifdef DEAL_II_WITH_T8CODE
#  include <deal.II/distributed/t8code_wrappers.h>
#endif // DEAL_II_WITH_T8CODE


namespace amr
{
#if defined(DEAL_II_WITH_T8CODE)
  using namespace dealii::internal::t8code;
#elif defined(DEAL_II_WITH_P4EST)
  using namespace dealii::internal::p4est;
#else
#  error DEAL_II_WITH_P4EST or DEAL_II_WITH_T8CODE required
#endif
} // namespace amr


#endif // dealii_p4est_wrappers_h
