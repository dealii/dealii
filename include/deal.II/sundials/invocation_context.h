// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_sundials_invocation_context_h
#define dealii_sundials_invocation_context_h

#include <deal.II/base/config.h>


#ifdef DEAL_II_WITH_SUNDIALS

#  include <sundials/sundials_types.h>

#  include <exception>


DEAL_II_NAMESPACE_OPEN


/**
 * A namespace for dealing with ODE solvers through the SUNDIALS package.
 */
namespace SUNDIALS
{
  namespace internal
  {
    /**
     * This struct contains references to the std::exception_ptr @p
     * pending_exception and SUNContext @p arkode_ctx for SUNDIALS 6.0 and
     * newer. This is required for situations when an object needs to pass
     * the context to the objects it uses.
     */
    struct InvocationContext
    {
      std::exception_ptr &pending_exception;
#  if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
      SUNContext &arkode_ctx;
#  endif
    };
  } // namespace internal
} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif


#endif
