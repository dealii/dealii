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
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
#ifndef MESQUITE_VECTOR3D_HPP
#define MESQUITE_VECTOR3D_HPP

#include "Mesquite.hpp"

#include <iosfwd>
#include <cassert>
#include <cstring>
#include <vector>

/*! \file Vector3D.hpp
  \brief Vector object with exactly 3 dimensions.

  This is as fast as a C array.

  \author Darryl Melander
  \author Thomas Leurent
*/
namespace MESQUITE_NS
{
   class Matrix3D;
   class MsqError;
   
   /*!
      \class Vector3D
      \brief Vector3D is the object that effeciently stores information about
      about three-deminsional vectors.  It is also the parent class of
      MsqVertex.      */
  class MESQUITE_EXPORT Vector3D
  {
  public:
    // Constructors
    Vector3D();
    Vector3D( double xyz );
    Vector3D( double x, double y, double z);
    Vector3D(const double xyz[3]);
    Vector3D(const Vector3D& to_copy);

    // *** virtual destructor *** Do not use for Vector3D, we need to keep
    // the deallocation of those objects very fast
    
    // Functions to get the coordinates
    double x() const;
    double y() const;
    double z() const;
    void get_coordinates(double& x, double& y, double& z) const;
    void get_coordinates(double xyz[3]) const;
    const double& operator[](size_t index) const; // 0-based
    
    // Functions to set the coordinates.
    void x(const double x);
    void y(const double y);
    void z(const double z);
    void set(const double x, const double y, const double z);
    void set(const double xyz[3]);
    void set(const Vector3D& to_copy);
    // Subscripts on non-consts both get and set coords
    double& operator[](size_t index); // 0-based
    Vector3D& operator=(const Vector3D &to_copy);
    Vector3D& operator=(const double &to_copy);
    
    // Functions to modify existing coordinates
    Vector3D operator-() const;  //- unary negation.
    Vector3D& operator*=(const double scalar);
    Vector3D& operator/=(const double scalar);
    Vector3D& operator*=(const Vector3D &rhs); //- cross product
    Vector3D& operator+=(const Vector3D &rhs);
    Vector3D& operator-=(const Vector3D &rhs);
    
    // Binary operators (like a+b).
    friend const Vector3D operator+(const Vector3D &lhs,
                              const Vector3D &rhs);
    friend const Vector3D operator-(const Vector3D &lhs,
                              const Vector3D &rhs);
    friend const Vector3D operator*(const Vector3D &lhs,
                              const double scalar); //!< lhs * scalar
    friend const Vector3D operator*(const double scalar,
                              const Vector3D &rhs); //!< scalar * rhs
    friend const Vector3D operator/(const Vector3D &lhs,
                              const double scalar); //- lhs / scalar
    friend double operator%(const Vector3D &v1,
                            const Vector3D &v2); //!< dot product
    friend double inner(const Vector3D v1[],
                        const Vector3D v2[], int n); //!< dot product for array
    friend double operator%(const double scalar,
			    const Vector3D &v2); //!< scalar * sum_i v2[i]
    friend double operator%(const Vector3D &v1,
                            const double scalar); //!< scalar * sum_i v1[i]
    friend const Vector3D operator*(const Vector3D &v1, 
                              const Vector3D &v2); //!< cross product
 
    //! \f$ v = A*x \f$
    friend void eqAx(Vector3D& v, const Matrix3D& A, const Vector3D& x);
    //! \f$ v += A*x \f$
    friend void plusEqAx(Vector3D& v, const Matrix3D& A, const Vector3D& x);
    //! \f$ v += A^T*x \f$
    friend void plusEqTransAx(Vector3D& v, const Matrix3D& A, const Vector3D& x);
    friend void eqTransAx(Vector3D& v, const Matrix3D& A, const Vector3D& x);
    
    // Comparison functions
    friend bool operator==(const Vector3D &lhs, const Vector3D &rhs);
    friend bool operator!=(const Vector3D &lhs, const Vector3D &rhs);
    static double distance_between(const Vector3D& p1,
                                   const Vector3D& p2);
    int within_tolerance_box(const Vector3D &compare_to,
                             double tolerance) const;
    //- Compare two Vector3Ds to see if they are spatially equal.  
    // Return TRUE if difference in x, y, and z are all within tolerance.
    // Essentially checks to see if 'this' lies within a box centered on
    // 'compare_to' with sides of length ('tolerance' * 2).

    // Length functions
    inline double length_squared() const;
    inline double length() const;
    friend  double length(const Vector3D* v,int n); //!< L2 norm for an array of Vector3Ds
    friend  double Linf(const Vector3D* v,int n); //!< L inf norm for array of Vector3Ds

    inline void set_length(const double new_length);
    inline void normalize();
    Vector3D operator~() const { return *this * (1.0/length()); }
    
    // Utility functions.  All angle functions work in radians.
    static double interior_angle(const Vector3D &a,
                                 const Vector3D &b,
                                 MsqError& err);
    //- Interior angle: acos((a%b)/(|a||b|))
    static Vector3D interpolate(const double param, const Vector3D &p1,
                                const Vector3D &p2);
    //- Interpolate between two points. Returns (1-param)*v1 + param*v2.

    const double* to_array() const 
      { return mCoords; }

    double* to_array() 
      { return mCoords; }

  protected:
    double mCoords[3];
  };

  // Constructors
  inline Vector3D::Vector3D() 
  {
    mCoords[0] = 0;
    mCoords[1] = 0;
    mCoords[2] = 0;
  }
  inline Vector3D::Vector3D( double x )
  {
    mCoords[0] = x;
    mCoords[1] = x;
    mCoords[2] = x;
  }
  inline Vector3D::Vector3D( double x,
                             double y,
                             double z) 
  {
    mCoords[0] = x;
    mCoords[1] = y;
    mCoords[2] = z;
  }
  inline Vector3D::Vector3D(const double xyz[3]) 
  { std::memcpy(mCoords, xyz, 3*sizeof(double)); }
  inline Vector3D::Vector3D(const Vector3D& to_copy) 
  { std::memcpy(mCoords, to_copy.mCoords, 3*sizeof(double)); }
  
  // Functions to get coordinates
  inline double Vector3D::x() const
  { return mCoords[0]; }
  inline double Vector3D::y() const
  { return mCoords[1]; }
  inline double Vector3D::z() const
  { return mCoords[2]; }
  inline void Vector3D::get_coordinates(double &x, double &y, double &z) const
  {
    x = mCoords[0]; 
    y = mCoords[1]; 
    z = mCoords[2];
  }
  inline void Vector3D::get_coordinates(double xyz[3]) const
  { std::memcpy(xyz, mCoords, 3*sizeof(double)); }
  inline const double& Vector3D::operator[](size_t index) const
  {
    return mCoords[index];
  }
  
  // Functions to set coordinates
  inline void Vector3D::x( const double x )
  { mCoords[0] = x; }
  inline void Vector3D::y( const double y )
  { mCoords[1] = y; }
  inline void Vector3D::z( const double z )
  { mCoords[2] = z; }
  inline void Vector3D::set(const double x,
                            const double y,
                            const double z)
  {
    mCoords[0] = x;
    mCoords[1] = y;
    mCoords[2] = z;
  }
  inline void Vector3D::set(const double xyz[3])
  { std::memcpy(mCoords, xyz, 3*sizeof(double)); }
  inline void Vector3D::set(const Vector3D& to_copy)
  { std::memcpy(mCoords, to_copy.mCoords, 3*sizeof(double)); }
  inline double& Vector3D::operator[](size_t index)
  { return mCoords[index]; }
   
  inline Vector3D& Vector3D::operator=(const Vector3D &to_copy)  
  {
    mCoords[0] = to_copy.mCoords[0];
    mCoords[1] = to_copy.mCoords[1];
    mCoords[2] = to_copy.mCoords[2];
//    memcpy(mCoords, to_copy.mCoords, 3*sizeof(double));
    return *this;
  }

  inline Vector3D& Vector3D::operator=(const double &to_copy)  
  {
    mCoords[0] = to_copy;
    mCoords[1] = to_copy;
    mCoords[2] = to_copy;
    return *this;
  }
  
  // Functions that modify existing coordinates
  inline Vector3D Vector3D::operator-() const
  {
    return Vector3D(-mCoords[0], -mCoords[1], -mCoords[2]);
  }
  inline Vector3D& Vector3D::operator*=(const double scalar)
  {
    mCoords[0] *= scalar;
    mCoords[1] *= scalar;
    mCoords[2] *= scalar;
    return *this;
  }
  //! divides each Vector3D entry by the given scalar.
  inline Vector3D& Vector3D::operator/=(const double scalar)
  {
    mCoords[0] /= scalar;
    mCoords[1] /= scalar;
    mCoords[2] /= scalar;
    return *this;
  }
  inline Vector3D& Vector3D::operator*=(const Vector3D &rhs)
  {
    double new_coords[3] = 
      {mCoords[1]*rhs.mCoords[2] - mCoords[2]*rhs.mCoords[1],
       mCoords[2]*rhs.mCoords[0] - mCoords[0]*rhs.mCoords[2],
       mCoords[0]*rhs.mCoords[1] - mCoords[1]*rhs.mCoords[0]
      };
    std::memcpy(mCoords, new_coords, 3*sizeof(double));
    return *this;
  }
  inline Vector3D& Vector3D::operator+=(const Vector3D &rhs)
  {
    mCoords[0] += rhs.mCoords[0];
    mCoords[1] += rhs.mCoords[1];
    mCoords[2] += rhs.mCoords[2];
    return *this;
  }
  inline Vector3D& Vector3D::operator-=(const Vector3D &rhs)
  {
    mCoords[0] -= rhs.mCoords[0];
    mCoords[1] -= rhs.mCoords[1];
    mCoords[2] -= rhs.mCoords[2];
    return *this;
  }

  // Binary operators
  inline const Vector3D operator+(const Vector3D &lhs,
                            const Vector3D &rhs)
  {
    return Vector3D(lhs.x() + rhs.x(),
                    lhs.y() + rhs.y(),
                    lhs.z() + rhs.z());
  }
  inline const Vector3D operator-(const Vector3D &lhs,
                            const Vector3D &rhs)
  {
    return Vector3D(lhs.x() - rhs.x(),
                    lhs.y() - rhs.y(),
                    lhs.z() - rhs.z());
  }
  inline const Vector3D operator*(const Vector3D &lhs,
                            const double scalar)
  {
    return Vector3D(lhs.x() * scalar,
                    lhs.y() * scalar,
                    lhs.z() * scalar);
  }
  inline const Vector3D operator*(const double scalar,
                            const Vector3D &rhs)
  {
    return Vector3D(rhs.x() * scalar,
                    rhs.y() * scalar,
                    rhs.z() * scalar);
  }
  inline const Vector3D operator/(const Vector3D &lhs,
                            const double scalar)
  {
    assert (scalar != 0);
    return Vector3D(lhs.x() / scalar,
                    lhs.y() / scalar,
                    lhs.z() / scalar);
  }
  inline double operator%(const Vector3D &lhs,
                          const Vector3D &rhs) // Dot Product
  {
    return( lhs.mCoords[0] * rhs.mCoords[0] +
            lhs.mCoords[1] * rhs.mCoords[1] +
            lhs.mCoords[2] * rhs.mCoords[2] );
  }

  /*! Dot product for arrays of Vector3Ds. see also operator% .*/ 
  inline double inner(const Vector3D lhs[],
                      const Vector3D rhs[], int n)
  {
    int i;
    double dot=0;
    for (i=0; i<n; ++i) 
      dot+= lhs[i].mCoords[0] * rhs[i].mCoords[0] +
            lhs[i].mCoords[1] * rhs[i].mCoords[1] +
            lhs[i].mCoords[2] * rhs[i].mCoords[2];
    return dot;
  }
  /*! Dot product for arrays of Vector3Ds. see also operator% .*/ 
  inline double inner(const std::vector<Vector3D>& lhs,
                      const std::vector<Vector3D>& rhs)
  {
    double dot = 0;
    assert(lhs.size() == rhs.size());
    for (size_t i = 0; i < lhs.size(); ++i)
      dot = lhs[i] % rhs[i];
    return dot;
  }

  inline double operator%(const double scalar,
                          const Vector3D &rhs) // Dot Product
  {
    return( scalar * (rhs.mCoords[0] + rhs.mCoords[1] + rhs.mCoords[2]) );
  }
  inline double operator%(const Vector3D &lhs,
                          const double scalar) // Dot Product
  {
    return( scalar * (lhs.mCoords[0] + lhs.mCoords[1] + lhs.mCoords[2]) );
  }
  inline const Vector3D operator*(const Vector3D &lhs, 
                            const Vector3D &rhs) // Cross Product
  {
    return Vector3D(lhs.mCoords[1]*rhs.mCoords[2]-lhs.mCoords[2]*rhs.mCoords[1],
                    lhs.mCoords[2]*rhs.mCoords[0]-lhs.mCoords[0]*rhs.mCoords[2],
                    lhs.mCoords[0]*rhs.mCoords[1]-lhs.mCoords[1]*rhs.mCoords[0]);
  }
  
  // output operator
  MESQUITE_EXPORT std::ostream& operator<<(std::ostream &s, const Mesquite::Vector3D &v);
  
  inline double Vector3D::distance_between(const Vector3D &p1,
                                           const Vector3D &p2)
  {
    Vector3D v = p2 - p1;
    return v.length();
  }
  inline int Vector3D::within_tolerance_box(const Vector3D &compare_to,
                                            double tolerance) const
  {
    return ((std::fabs(this->mCoords[0] - compare_to.mCoords[0]) < tolerance) &&
            (std::fabs(this->mCoords[1] - compare_to.mCoords[1]) < tolerance) &&
            (std::fabs(this->mCoords[2] - compare_to.mCoords[2]) < tolerance));
  }
  
  // Length functions
  inline double Vector3D::length_squared() const
  {
    return (mCoords[0]*mCoords[0] +
            mCoords[1]*mCoords[1] +
            mCoords[2]*mCoords[2]);
  }
  inline double Vector3D::length() const
  {
    return std::sqrt(mCoords[0]*mCoords[0] +
                          mCoords[1]*mCoords[1] +
                          mCoords[2]*mCoords[2]);
  }
  
  inline double inner_product( const Vector3D* v1, const Vector3D* v2, size_t n )
  {
    double result = 0.0;
    const Vector3D* const end = v1 + n;
    while (v1 < end) {
      result += *v1 % *v2;
      ++v1;
      ++v2;
    }
    return result;
  }

  inline double length_squared( const Vector3D*  v, int n )
  {
    double sum = 0.0;
    for (int i = 0; i < n; ++i)
      sum += v[i].length_squared();
    return sum;
  }
  inline double length_squared( const std::vector<Vector3D>& v )
  {
    double sum = 0.0;
    for (size_t i = 0; i < v.size(); ++i)
      sum += v[i].length_squared();
    return sum;
  }

  inline double length(const Vector3D*  v,int n) // norm for an array of Vector3Ds
  {
    return std::sqrt( length_squared( v, n ) );
  }
  inline double length( const std::vector<Vector3D>& v )
  {
    return std::sqrt( length_squared( v ) );
  }

  inline double Linf(const Vector3D*  v,int n) // max entry for an array of Vector3Ds
  {
    double max=0;  
    //loop over the length of the array
    for(int i=0;i<n;++i){
      if ( max < std::fabs(v[i][0]) )   max=std::fabs(v[i][0]) ;
      if ( max < std::fabs(v[i][1]) )   max=std::fabs(v[i][1]) ;
      if ( max < std::fabs(v[i][2]) )   max=std::fabs(v[i][2]) ;
    }
    //return the value of the largest entry in the array
    return max;
  }

  inline double Linf( const std::vector<Vector3D>& v ) // max entry for an array of Vector3Ds
  {
    double max=0;  
    //loop over the length of the array
    for(size_t i=0;i<v.size();++i){
      if ( max < std::fabs(v[i][0]) )   max=std::fabs(v[i][0]) ;
      if ( max < std::fabs(v[i][1]) )   max=std::fabs(v[i][1]) ;
      if ( max < std::fabs(v[i][2]) )   max=std::fabs(v[i][2]) ;
    }
    //return the value of the largest entry in the array
    return max;
  }

  inline void Vector3D::set_length(const double new_length)
  {
    double factor = new_length / length();
    *this *= factor;
  }
  inline void Vector3D::normalize()
  { set_length(1.0); }

  // Utility functions.
  inline Vector3D Vector3D::interpolate(const double param,
                                        const Vector3D &p1,
                                        const Vector3D &p2)
  {
    return (1-param)*p1 + param*p2;
  }

  inline bool operator==( const Vector3D& v1, const Vector3D&v2 )
    { return v1.mCoords[0] == v2.mCoords[0] &&
             v1.mCoords[1] == v2.mCoords[1] &&
             v1.mCoords[2] == v2.mCoords[2]; }

  inline bool operator!=( const Vector3D& v1, const Vector3D&v2 )
    { return v1.mCoords[0] != v2.mCoords[0] ||
             v1.mCoords[1] != v2.mCoords[1] ||
             v1.mCoords[2] != v2.mCoords[2]; }

} // namespace Mesquite

#endif
