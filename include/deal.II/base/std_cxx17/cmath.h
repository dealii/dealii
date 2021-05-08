// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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



  inline long double
  cyl_bessel_jl(long double nu, long double x)
  {
    return boost::math::cyl_bessel_j(nu, x);
  }

#else
  using std::cyl_bessel_j;
  using std::cyl_bessel_jf;
  using std::cyl_bessel_jl;
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



  inline long double
  legendre(unsigned int l, long double x)
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



  inline long double
  legendrel(unsigned int l, long double x)
  {
    Assert(static_cast<int>(l) >= 0,
           ExcIndexRange(l, 0, std::numeric_limits<int>::max()));
    return boost::math::legendre_p(static_cast<int>(l), x);
  }

#else
  using std::legendre;
  using std::legendref;
  using std::legendrel;
#endif
} // namespace std_cxx17


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_cxx17_cmath_h
