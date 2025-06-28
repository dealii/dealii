// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_differentiation_sd_h
#define dealii_differentiation_sd_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/differentiation/sd/symengine_math.h>
#  include <deal.II/differentiation/sd/symengine_number_traits.h>
#  include <deal.II/differentiation/sd/symengine_number_types.h>
#  include <deal.II/differentiation/sd/symengine_optimizer.h>
#  include <deal.II/differentiation/sd/symengine_product_types.h>
#  include <deal.II/differentiation/sd/symengine_scalar_operations.h>
#  include <deal.II/differentiation/sd/symengine_tensor_operations.h>
#  include <deal.II/differentiation/sd/symengine_types.h>
#  include <deal.II/differentiation/sd/symengine_utilities.h>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  /**
   * Wrappers for symbolic differentiation libraries. Currently there is support
   * for the following libraries:
   *   - SymEngine
   *
   * @ingroup auto_symb_diff
   */
  namespace SD
  {
    /**
     * This namespace defines the classes and functions that help provide a
     * structured interface to symbolic numbers and operations.
     */
    namespace internal
    {}
  } // namespace SD
} // namespace Differentiation

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif // dealii_differentiation_sd_h
