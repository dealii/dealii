// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_differentiation_ad_adolc_math_h
#define dealii_differentiation_ad_adolc_math_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_ADOLC

#  include <adolc/adouble.h> // Taped double
#  include <adolc/adtl.h>    // Tapeless double
#  include <adolc/internal/adolc_settings.h>
#  include <adolc/internal/adubfunc.h> // Taped double math functions

#  include <cmath>


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

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_ADOLC

#endif
