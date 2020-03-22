/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov, 
    kraftche@cae.wisc.edu     
   
  ***************************************************************** */
/*!
  \file   Exponent.cpp
  \brief 

  \author Jason Kraftcheck
  \date   2005-5-2
*/

#include "Exponent.hpp"

namespace MESQUITE_NS {



Exponent::constMemberPtr Exponent::get_func_ptr( double exponent )
{
  if (exponent == 0.0)
    return &Exponent::pow0;
  else if (exponent == 1.0)
    return &Exponent::pow1;
  else if (exponent == 0.5)
    return &Exponent::squareRoot;
  else if (exponent == 1./3.)
    return &Exponent::cubeRoot;
  else if (exponent == -1./3.)
    return &Exponent::invCubeRoot;
  else if (exponent == 2./3.)
    return &Exponent::powTwoThirds;
  else if (exponent == -2./3.)
    return &Exponent::invTwoThirds;
  else if (exponent == 2.0)
    return &Exponent::pow2;
  else if (exponent == -1.0)
    return &Exponent::inverse;
  else if (exponent == -0.5)
    return &Exponent::invSquareRoot;
  else if (exponent == 1.5)
    return &Exponent::powThreeHalves;
  else if (exponent == -2.0)
    return &Exponent::invSquare;
  else if (std::floor(exponent) == exponent)
  {
    if (exponent > 0.0)
      return &Exponent::powPositiveInt;
    else
      return &Exponent::powNegativeInt;
  }
  else 
    return &Exponent::std_pow;
}

void Exponent::set_exponent( double exponent )
{
  mExponent = exponent;
  funcPointer = get_func_ptr( exponent );
}

double Exponent::pow0( double   ) const          { return 1.0;   }
double Exponent::pow1( double x ) const          { return x;     }
double Exponent::pow2( double x ) const          { return x * x; }
double Exponent::squareRoot( double x ) const    { return std::sqrt( x ); }
double Exponent::cubeRoot( double x ) const      { return Mesquite::cbrt( x ); }
double Exponent::invCubeRoot( double x ) const   { return 1.0/Mesquite::cbrt( x ); }
double Exponent::powTwoThirds( double x ) const  { return Mesquite::cbrt_sqr(x); }
double Exponent::invTwoThirds( double x ) const  { return 1.0 / Mesquite::cbrt_sqr(x); }
double Exponent::std_pow( double x ) const       { return std::pow( x, mExponent ); }
double Exponent::inverse( double x ) const       { return 1.0 / x; }
double Exponent::invSquareRoot( double x ) const { return 1.0 / std::sqrt(x); }
double Exponent::powThreeHalves( double x ) const{ return x*x*x / std::sqrt(x); }
double Exponent::invSquare( double x ) const     { return 1.0 / (x*x); }

double Exponent::powPositiveInt( double x ) const
{
  double result = x;
  for (int i = (int)mExponent - 1; i > 0; --i)
    result *= x;
  return result;
}

double Exponent::powNegativeInt( double x ) const
{
  double result = x;
  for (int i = (-(int)mExponent) - 1; i > 0; --i)
    result *= x;
  return 1.0/result;
}


} // namespace Mesquite
