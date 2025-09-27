// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_differentiation_ad_sacado_math_h
#define dealii_differentiation_ad_sacado_math_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_SACADO

#  include <deal.II/base/numbers.h>

#  include <deal.II/differentiation/ad/ad_number_traits.h>
#  include <deal.II/differentiation/ad/sacado_number_types.h>


#  ifndef DOXYGEN

DEAL_II_NAMESPACE_OPEN
/**
 * Implementation of the error function for real-valued Sacado numbers.
 */
template <
  typename ADNumberType,
  typename = std::enable_if_t<
    dealii::Differentiation::AD::is_sacado_number<ADNumberType>::value &&
    dealii::Differentiation::AD::is_real_valued_ad_number<ADNumberType>::value>>
inline ADNumberType
erf(ADNumberType x)
{
  // Reference:
  // Handbook of Mathematical Functions: with Formulas, Graphs, and
  // Mathematical Tables Abramowitz, M. and Stegun, I. Dover Books on
  // Mathematics. 1972.
  //
  // Current implementation source:
  // https://www.johndcook.com/blog/cpp_erf/
  // https://www.johndcook.com/blog/2009/01/19/stand-alone-error-function-erf/
  //
  // Note: This implementation has a reported maximum round-off error
  // of 1.5e-7.

  // Constants
  const double a1 = 0.254829592;
  const double a2 = -0.284496736;
  const double a3 = 1.421413741;
  const double a4 = -1.453152027;
  const double a5 = 1.061405429;
  const double p  = 0.3275911;

  // Save the sign of x
  const bool neg_val =
    (x < dealii::internal::NumberType<ADNumberType>::value(0.0) ? true : false);
  x = std::fabs(x);

  // Abramowitz1972a equation 7.1.26
  const ADNumberType t = 1.0 / (1.0 + p * x);
  const ADNumberType y =
    1.0 -
    (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);

  if (!neg_val)
    return y;
  else
    return -y;
}


/**
 * Implementation of the complementary error function for Sacado numbers.
 */
template <typename ADNumberType,
          typename = std::enable_if_t<
            dealii::Differentiation::AD::is_sacado_number<ADNumberType>::value>>
inline ADNumberType
erfc(const ADNumberType &x)
{
  return 1.0 - std::erf(x);
}

DEAL_II_NAMESPACE_CLOSE

#  endif // DOXYGEN

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_SACADO

#endif
