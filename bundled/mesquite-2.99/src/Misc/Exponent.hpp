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
  \file   Exponent.hpp
  \brief 

  \author Jason Kraftcheck
  \date   2005-5-2
*/

#ifndef MESQUITE_POWER_HPP
#define MESQUITE_POWER_HPP

#include "Mesquite.hpp"

namespace MESQUITE_NS {

class Exponent
{
  public:
    
    typedef double (Exponent::*constMemberPtr)(double) const;
    static constMemberPtr get_func_ptr( double exponent );
  
    Exponent( ) : funcPointer( 0 ) 
      {}
  
    explicit Exponent( double exponent ) 
      : mExponent(exponent), 
        funcPointer( get_func_ptr( exponent ) )
      {}
    
    inline double raise( double value ) const
      { return (this->*funcPointer)( value ); }
    
    void set_exponent( double exponent );
    
    inline Exponent& operator=( double d )
      { set_exponent(d); return *this; }
    
    //inline operator double () const
    //  { return mExponent; }
    inline double value() const
      { return mExponent; }
    
    double pow0( double x ) const;
    double pow1( double x ) const;
    double squareRoot( double x ) const;
    double cubeRoot( double x ) const;
    double invCubeRoot( double x ) const;
    double powTwoThirds( double x ) const;
    double invTwoThirds( double x ) const;
    double pow2( double x ) const;
    double powPositiveInt( double x ) const;
    double std_pow( double x ) const;
    double inverse( double x ) const;
    double invSquareRoot( double x ) const;
    double powThreeHalves( double x ) const;
    double invSquare( double x ) const;
    double powNegativeInt( double x ) const;

  private:
  
    double mExponent;
    constMemberPtr funcPointer;
};

} // namespace Mesquite

#endif
