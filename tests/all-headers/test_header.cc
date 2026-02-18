// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// include header twice to see if the include guards are set correctly
#include HEADER
#include HEADER

#if !defined(DEAL_II_NAMESPACE_OPEN) && !defined(dealii_revision_h)
#  error "HEADER does not include config.h."
#endif


int
main()
{
  return 0;
}
