// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_sundials_types_h
#define dealii_sundials_types_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS
#  include <sundials/sundials_types.h>

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
/**
 * Alias for the real type used by SUNDIALS.
 */
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
  using realtype = ::sunrealtype;
#  else
  using realtype = ::realtype;
#  endif
} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SUNDIALS
#endif // dealii_sundials_types_h
