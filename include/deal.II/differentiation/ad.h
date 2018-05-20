// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii_differentiation_ad_h
#define dealii_differentiation_ad_h

#include <deal.II/base/config.h>

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
   * Wrappers for automatic differentiation libraries. Currently there is support
   * for the following libraries:
   *   - Adol-C
   *   - Sacado (a component of Trilinos)
   *
   * @ingroup auto_symb_diff
   */
  namespace AD
  {}
} // namespace Differentiation

DEAL_II_NAMESPACE_CLOSE

#endif
