// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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

#ifdef DEAL_II_HAVE_CXX17_BESSEL_FUNCTIONS
#  include <cmath>
#else
#  include <boost/math/special_functions/bessel.hpp>
#endif


DEAL_II_NAMESPACE_OPEN
namespace std_cxx17
{
#ifndef DEAL_II_HAVE_CXX17_BESSEL_FUNCTIONS

  inline double
  cyl_bessel_j(double x, double y)
  {
    return boost::math::cyl_bessel_j(x, y);
  }



  inline float
  cyl_bessel_jf(float x, float y)
  {
    return boost::math::cyl_bessel_j(x, y);
  }



  inline long double
  cyl_bessel_jl(long double x, long double y)
  {
    return boost::math::cyl_bessel_j(x, y);
  }

#else
  using std::cyl_bessel_j;
  using std::cyl_bessel_jf;
  using std::cyl_bessel_jl;
#endif
} // namespace std_cxx17
DEAL_II_NAMESPACE_CLOSE

#endif // dealii_cxx17_cmath_h
