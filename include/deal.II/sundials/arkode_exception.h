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


#ifndef dealii_sundials_arkode_exception_h
#define dealii_sundials_arkode_exception_h

#include <deal.II/base/exception_macros.h>
#include <deal.II/base/exceptions.h>


#ifdef DEAL_II_WITH_SUNDIALS


DEAL_II_NAMESPACE_OPEN


/**
 * A namespace for dealing with ODE solvers through the SUNDIALS package.
 */
namespace SUNDIALS
{
  /**
   * Handle ARKode exceptions.
   */
  DeclException1(ExcARKodeError,
                 int,
                 << "One of the SUNDIALS ARKode internal functions "
                 << " returned a negative error code: " << arg1
                 << ". Please consult SUNDIALS manual.");
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
