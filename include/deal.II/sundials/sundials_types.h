// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_sundials_types_h
#define dealii_sundials_types_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS
#  include <sundials/sundials_types.h>

#endif // DEAL_II_WITH_SUNDIALS

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_SUNDIALS
namespace SUNDIALS
{
/**
 * Alias for the bool and real types used by SUNDIALS.
 */
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
  using booltype = sunbooleantype;
  using realtype = ::sunrealtype;
#  else
  using booltype = booleantype;
  using realtype = ::realtype;
#  endif
} // namespace SUNDIALS

#endif // DEAL_II_WITH_SUNDIALS

DEAL_II_NAMESPACE_CLOSE
#endif // dealii_sundials_types_h
