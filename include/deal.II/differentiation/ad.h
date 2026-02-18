// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_differentiation_ad_h
#define dealii_differentiation_ad_h

#include <deal.II/base/config.h>

#include <deal.II/differentiation/ad/ad_helpers.h>
#include <deal.II/differentiation/ad/ad_number_traits.h>
#include <deal.II/differentiation/ad/ad_number_types.h>
#include <deal.II/differentiation/ad/adolc_math.h>
#include <deal.II/differentiation/ad/adolc_number_types.h>
#include <deal.II/differentiation/ad/adolc_product_types.h>
#include <deal.II/differentiation/ad/sacado_math.h>
#include <deal.II/differentiation/ad/sacado_number_types.h>
#include <deal.II/differentiation/ad/sacado_product_types.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace that encapsulates various classes and helper functions related
 * to automatic and symbolic differentiation.
 *
 * @ingroup auto_symb_diff
 */
namespace Differentiation
{
  /**
   * Wrappers for automatic differentiation libraries. Currently there is
   * support for the following libraries:
   *   - ADOL-C
   *   - Sacado (a component of Trilinos)
   *
   * @ingroup auto_symb_diff
   */
  namespace AD
  {}
} // namespace Differentiation

DEAL_II_NAMESPACE_CLOSE

#endif
