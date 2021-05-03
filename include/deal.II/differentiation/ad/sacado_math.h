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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_differentiation_ad_sacado_math_h
#define dealii_differentiation_ad_sacado_math_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_SACADO

#  include <deal.II/base/numbers.h>

#  include <deal.II/differentiation/ad/ad_number_traits.h>
#  include <deal.II/differentiation/ad/sacado_number_types.h>


#  ifndef DOXYGEN

/**
 * Define some missing fundamental math functions
 */
namespace std
{
  /**
   * Implementation of the error function for real-valued Sacado numbers.
   */
  template <
    typename ADNumberType,
    typename = typename std::enable_if<
      dealii::Differentiation::AD::is_sacado_number<ADNumberType>::value &&
      dealii::Differentiation::AD::is_real_valued_ad_number<
        ADNumberType>::value>::type>
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
      (x < dealii::internal::NumberType<ADNumberType>::value(0.0) ? true :
                                                                    false);
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
  template <
    typename ADNumberType,
    typename = typename std::enable_if<
      dealii::Differentiation::AD::is_sacado_number<ADNumberType>::value>::type>
  inline ADNumberType
  erfc(const ADNumberType &x)
  {
    return 1.0 - std::erf(x);
  }

} // namespace std

#  endif // DOXYGEN

#endif // DEAL_II_TRILINOS_WITH_SACADO

#endif
