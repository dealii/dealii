// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/symbolic_function.templates.h>

DEAL_II_NAMESPACE_OPEN
#ifdef DEAL_II_WITH_SYMENGINE
namespace Functions
{
  // explicit instantiations
#  include "base/symbolic_function.inst"
} // namespace Functions
#endif
DEAL_II_NAMESPACE_CLOSE
