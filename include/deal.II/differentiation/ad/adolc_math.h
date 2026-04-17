// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_differentiation_ad_adolc_math_h
#define dealii_differentiation_ad_adolc_math_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_ADOLC

#  include <adolc/adouble.h> // Taped double
#  include <adolc/adtl.h>    // Tapeless double
#  include <adolc/internal/adolc_settings.h>
#  include <adolc/internal/adubfunc.h> // Taped double math functions

#  include <cmath>

#endif // DEAL_II_WITH_ADOLC

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_ADOLC
#  ifndef DOXYGEN
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

#  endif // DOXYGEN

#endif // DEAL_II_WITH_ADOLC

DEAL_II_NAMESPACE_CLOSE

#endif
