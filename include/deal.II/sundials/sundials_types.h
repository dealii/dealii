// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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
