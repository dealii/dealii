// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
#ifndef dealii_cxx17_cmath_h
#define dealii_cxx17_cmath_h

#include <deal.II/base/config.h>

#if defined(DEAL_II_HAVE_CXX17_BESSEL_FUNCTIONS) || \
  defined(DEAL_II_HAVE_CXX17_LEGENDRE_FUNCTIONS)
#  include <cmath>
#endif

#ifndef DEAL_II_HAVE_CXX17_BESSEL_FUNCTIONS
#  include <boost/math/special_functions/bessel.hpp>
#endif

#ifndef DEAL_II_HAVE_CXX17_LEGENDRE_FUNCTIONS
#  include <deal.II/base/exceptions.h>

#  include <boost/math/special_functions/legendre.hpp>

#  include <limits>
#endif


DEAL_II_NAMESPACE_OPEN

namespace std_cxx17
{
#ifndef DEAL_II_HAVE_CXX17_BESSEL_FUNCTIONS

  inline double
  cyl_bessel_j(double nu, double x)
  {
    return boost::math::cyl_bessel_j(nu, x);
  }



  inline float
  cyl_bessel_jf(float nu, float x)
  {
    return boost::math::cyl_bessel_j(nu, x);
  }

#else
  using std::cyl_bessel_j;
  using std::cyl_bessel_jf;
#endif

#ifndef DEAL_II_HAVE_CXX17_LEGENDRE_FUNCTIONS

  inline double
  legendre(unsigned int l, double x)
  {
    Assert(static_cast<int>(l) >= 0,
           ExcIndexRange(l, 0, std::numeric_limits<int>::max()));
    return boost::math::legendre_p(static_cast<int>(l), x);
  }



  inline float
  legendre(unsigned int l, float x)
  {
    Assert(static_cast<int>(l) >= 0,
           ExcIndexRange(l, 0, std::numeric_limits<int>::max()));
    return boost::math::legendre_p(static_cast<int>(l), x);
  }



  inline float
  legendref(unsigned int l, float x)
  {
    Assert(static_cast<int>(l) >= 0,
           ExcIndexRange(l, 0, std::numeric_limits<int>::max()));
    return boost::math::legendre_p(static_cast<int>(l), x);
  }

#else
  using std::legendre;
  using std::legendref;
#endif
} // namespace std_cxx17


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_cxx17_cmath_h
