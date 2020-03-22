/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file SymMatrix3D.hpp
 *  \brief Symetric 3x3 Matrix
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_SYM_MATRIX_3D_HPP
#define MSQ_SYM_MATRIX_3D_HPP

#include "Mesquite.hpp"
#include "Vector3D.hpp"

namespace MESQUITE_NS {

class MESQUITE_EXPORT SymMatrix3D {
private:
  double d_[6];

public:
  enum Term {
    T00 = 0,
    T01 = 1,
    T02 = 2,
    T10 = T01,
    T11 = 3,
    T12 = 4,
    T20 = T02,
    T21 = T12,
    T22 = 5
  };
  
  inline static Term term( unsigned r, unsigned c )
    { return (Term)(r <= c ? 3*r - r*(r+1)/2 + c : 3*c - c*(c+1)/2 + r); }
  
  SymMatrix3D() {}
  
  SymMatrix3D( double diagonal_value ) {
    d_[T00] = d_[T11] = d_[T22] = diagonal_value;
    d_[T01] = d_[T02] = d_[T12] = 0.0;
  }
  
  SymMatrix3D( double t00, double t01, double t02,
                           double t11, double t12,
                                       double t22 ) 
  {
    d_[T00] = t00;
    d_[T01] = t01;
    d_[T02] = t02;
    d_[T11] = t11;
    d_[T12] = t12;
    d_[T22] = t22;
  }
    
    /**\brief Outer product */
  SymMatrix3D( const Vector3D& u )
  {
    d_[T00] = u[0] * u[0];
    d_[T01] = u[0] * u[1];
    d_[T02] = u[0] * u[2];
    d_[T11] = u[1] * u[1];
    d_[T12] = u[1] * u[2];
    d_[T22] = u[2] * u[2];
  }
    
 
  double& operator[]( unsigned t )       { return d_[t]; }
  double  operator[]( unsigned t ) const { return d_[t]; }
  
  double& operator()( unsigned short r, unsigned short c )
    { return d_[term(r,c)]; }
  double  operator()( unsigned short r, unsigned short c ) const
    { return d_[term(r,c)]; }
  
  inline SymMatrix3D& operator+=( const SymMatrix3D& other );
  inline SymMatrix3D& operator-=( const SymMatrix3D& other );
  inline SymMatrix3D& operator*=( double scalar );
  inline SymMatrix3D& operator/=( double scalar );
};

inline SymMatrix3D operator-( const SymMatrix3D& m )
{
  return SymMatrix3D( 
    -m[SymMatrix3D::T00], -m[SymMatrix3D::T01], -m[SymMatrix3D::T02],
                          -m[SymMatrix3D::T11], -m[SymMatrix3D::T12],
                                                -m[SymMatrix3D::T22] );
}

inline SymMatrix3D& SymMatrix3D::operator+=( const SymMatrix3D& other )
{
  d_[0] += other.d_[0];
  d_[1] += other.d_[1];
  d_[2] += other.d_[2];
  d_[3] += other.d_[3];
  d_[4] += other.d_[4];
  d_[5] += other.d_[5];
  return *this;
}

inline SymMatrix3D& SymMatrix3D::operator-=( const SymMatrix3D& other )
{
  d_[0] -= other.d_[0];
  d_[1] -= other.d_[1];
  d_[2] -= other.d_[2];
  d_[3] -= other.d_[3];
  d_[4] -= other.d_[4];
  d_[5] -= other.d_[5];
  return *this;
}

inline SymMatrix3D& SymMatrix3D::operator*=( double s )
{
  d_[0] *= s;
  d_[1] *= s;
  d_[2] *= s;
  d_[3] *= s;
  d_[4] *= s;
  d_[5] *= s;
  return *this;
}

inline SymMatrix3D& SymMatrix3D::operator/=( double s )
{
  d_[0] /= s;
  d_[1] /= s;
  d_[2] /= s;
  d_[3] /= s;
  d_[4] /= s;
  d_[5] /= s;
  return *this;
}


inline SymMatrix3D operator+( const SymMatrix3D& a, const SymMatrix3D& b )
  { SymMatrix3D r(a); r += b; return r; }
inline SymMatrix3D operator-( const SymMatrix3D& a, const SymMatrix3D& b )
  { SymMatrix3D r(a); r -= b; return r; }
inline SymMatrix3D operator*( const SymMatrix3D& a, double s )
  { SymMatrix3D r(a); r *= s; return r; }
inline SymMatrix3D operator*( double s, const SymMatrix3D& a )
  { SymMatrix3D r(a); r *= s; return r; }
inline SymMatrix3D operator/( const SymMatrix3D& a, double s )
  { SymMatrix3D r(a); r /= s; return r; }
inline SymMatrix3D operator/( double s, const SymMatrix3D& a )
  { SymMatrix3D r(a); r /= s; return r; }

inline Vector3D operator*( const Vector3D& v, const SymMatrix3D& m )
{
  return Vector3D( v[0]*m[0] + v[1]*m[1] + v[2]*m[2],
                   v[0]*m[1] + v[1]*m[3] + v[2]*m[4],
                   v[0]*m[2] + v[1]*m[4] + v[2]*m[5] );
}
inline Vector3D operator*( const SymMatrix3D& m, const Vector3D& v )
{
  return v * m;
}

/** Calculate the outer product of a vector with itself */
inline SymMatrix3D outer( const Vector3D& v )
{
  return SymMatrix3D( v[0]*v[0], v[0]*v[1], v[0]*v[2],
                                 v[1]*v[1], v[1]*v[2],
                                            v[2]*v[2] );
}

/** Given to vectors u and v, calculate the symmetric matrix
 *  equal to outer(u,v) + transpose(outer(u,v))
 *  equal to outer(v,u) + transpose(outer(v,u))
 */
inline SymMatrix3D outer_plus_transpose( const Vector3D& u,
                                         const Vector3D& v )
{
  return SymMatrix3D( 2*u[0]*v[0], u[0]*v[1] + u[1]*v[0], u[0]*v[2] + u[2]*v[0],
                                   2*u[1]*v[1]          , u[1]*v[2] + u[2]*v[1],
                                                          2*u[2]*v[2] );
}

inline const SymMatrix3D& transpose( const SymMatrix3D& a )
  { return a; }

inline double det( const SymMatrix3D& a )
{
  return a[0]*a[3]*a[5] + 2.0*a[1]*a[2]*a[4]
       - a[0]*a[4]*a[4] - a[3]*a[2]*a[2] - a[5]*a[1]*a[1];
}

inline SymMatrix3D inverse( const SymMatrix3D& a )
{
  SymMatrix3D result( a[3]*a[5] - a[4]*a[4],
                      a[2]*a[4] - a[1]*a[5],
                      a[1]*a[4] - a[2]*a[3],
                      a[0]*a[5] - a[2]*a[2],
                      a[1]*a[2] - a[0]*a[4],
                      a[0]*a[3] - a[1]*a[1] );
  result /= det( a );
  return result;
}

inline double Frobenius_2( const SymMatrix3D& a )
{
  return   a[0]*a[0] +
         2*a[1]*a[1] +
         2*a[2]*a[2] +
           a[3]*a[3] +
         2*a[4]*a[5] +
           a[5]*a[5];
}

inline double Frobenius( const SymMatrix3D& a )
  { return std::sqrt( Frobenius_2(a) ); }


} // namespace Mesquite

#endif
