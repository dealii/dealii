// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2021 by the deal.II authors
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

#ifndef dealii_differentiation_ad_adolc_math_h
#define dealii_differentiation_ad_adolc_math_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_ADOLC

#  include <adolc/adouble.h> // Taped double
#  include <adolc/adtl.h>    // Tapeless double
#  include <adolc/internal/adolc_settings.h>
#  include <adolc/internal/adubfunc.h> // Taped double math functions


#  ifndef DOXYGEN

DEAL_II_NAMESPACE_OPEN
/**
 * Implementation of the complementary error function for adol-c adouble
 * numbers.
 */
#    ifdef DEAL_II_ADOLC_WITH_ATRIG_ERF
inline adouble
erfc(const adouble &x)
{
  // Make things work with AD types
  using std::erf;
  return 1.0 - erf(x);
}

/**
 * Implementation of the complementary error function for adol-c adtl::adouble
 * numbers.
 */
inline adtl::adouble
erfc(const adtl::adouble &x)
{
  // Make things work with AD types
  using std::erf;
  return 1.0 - erf(x);
}
#    endif


/**
 * For improved genericity, implement abs() in terms of fabs() for adouble.
 */
inline adouble
abs(const adouble &x)
{
  return fabs(static_cast<const badouble &>(x));
}

/**
 * Same idea: implement abs() in terms of adtl::fabs().
 */
inline adtl::adouble
abs(const adtl::adouble &x)
{
  return adtl::fabs(x);
}

DEAL_II_NAMESPACE_CLOSE
#  endif // DOXYGEN

#endif // DEAL_II_WITH_ADOLC

#endif
