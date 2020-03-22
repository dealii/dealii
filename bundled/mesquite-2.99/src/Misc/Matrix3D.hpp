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
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 18-Dec-02 at 11:08:22
//  LAST-MOD: 27-May-04 at 14:48:56 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file Matrix3D.hpp

3*3 Matric class, row-oriented, 0-based [i][j] indexing.

 \author Thomas Leurent
 
*/
// DESCRIP-END.
//



#ifndef Matrix3D_hpp
#define Matrix3D_hpp

#include <iostream>
#include <sstream>
#include <cstdlib>

#include "Mesquite.hpp"
#include "Vector3D.hpp"
#include "SymMatrix3D.hpp"

namespace MESQUITE_NS
{

  /*! \class Matrix3D
      \brief 3*3 Matric class, row-oriented, 0-based [i][j] indexing.

      Since the size of the object is fixed at compile time, the Matrix3D
      object is as fast as a double[9] array.
  */
  class MESQUITE_EXPORT Matrix3D 
  {
  protected:
    double v_[9];                  

   
    void copy(const double*  v)
    { 
      v_[0] = v[0];
      v_[1] = v[1];
      v_[2] = v[2];
      v_[3] = v[3];
      v_[4] = v[4];
      v_[5] = v[5];
      v_[6] = v[6];
      v_[7] = v[7];
      v_[8] = v[8];
    }

    void set(double val)
    {
      v_[0]=val;  v_[1]=val;  v_[2]=val;
      v_[3]=val;  v_[4]=val;  v_[5]=val;
      v_[6]=val;  v_[7]=val;  v_[8]=val;
    }

    inline void set_values(const char *s);
    
  public:

    // constructors
    //! Default constructor sets all entries to 0. 
    Matrix3D()
    {
      zero();
    }
    
    Matrix3D(const Matrix3D &A)
    {
      copy(A.v_);
    }

    //! sets all entries of the matrix to value.
    Matrix3D(double value)
    {
      set(value);
    }
    
    Matrix3D( double a00, double a01, double a02,
              double a10, double a11, double a12,
              double a20, double a21, double a22 )
    {
      v_[0] = a00; v_[1] = a01; v_[2] = a02;
      v_[3] = a10; v_[4] = a11; v_[5] = a12;
      v_[6] = a20; v_[7] = a21; v_[8] = a22;
    }
    
    Matrix3D( const Vector3D& col1,
              const Vector3D& col2,
              const Vector3D& col3 )
    {
      set_column( 0, col1 );
      set_column( 1, col2 );
      set_column( 2, col3 );
    }
    
    Matrix3D( double radians, const Vector3D& axis )
    {
      Vector3D v(axis);
      v.normalize();
      const double c = std::cos( radians );
      const double s = std::sin( radians );
      v_[0] =  c      + (1.0 - c) * v[0]*v[0];
      v_[1] = -v[2]*s + (1.0 - c) * v[0]*v[1];
      v_[2] =  v[1]*s + (1.0 - c) * v[0]*v[2];
      v_[3] =  v[2]*s + (1.0 - c) * v[0]*v[1];
      v_[4] =  c      + (1.0 - c) * v[1]*v[1];
      v_[5] = -v[0]*s + (1.0 - c) * v[1]*v[2];
      v_[6] = -v[1]*s + (1.0 - c) * v[0]*v[2];
      v_[7] =  v[0]*s + (1.0 - c) * v[1]*v[2];
      v_[8] =  c      + (1.0 - c) * v[2]*v[2];
    }
      

    //! sets matrix entries to values in array.
    //! \param v is an array of 9 doubles. 
    Matrix3D(const double* v)
    {
      copy(v);
    }

    //! for test purposes, matrices can be instantiated as
    //! \code Matrix3D A("3 2 1  4 5 6  9 8 7"); \endcode
    Matrix3D(const char *s)
    {
      set_values(s);
    }
    
    Matrix3D( const SymMatrix3D& m )
    {
      *this = m;
    }

    // destructor
    ~Matrix3D() { }

    // assignments
    Matrix3D& operator=(const Matrix3D &A)
    {
      copy(A.v_);
      return *this;
    }
    
    Matrix3D& operator=( const SymMatrix3D& m )
    {
      v_[0]         = m[0];
      v_[1] = v_[3] = m[1];
      v_[2] = v_[6] = m[2];
      v_[4]         = m[3];
      v_[5] = v_[7] = m[4];
      v_[8]         = m[5];
      return *this;
    }
    
        
    Matrix3D& operator=(double scalar)
    { 
      set(scalar); 
      return *this;
    }

    //! for test purposes, matrices can be assigned as follows
    //! \code A = "3 2 1  4 5 6  9 8 7"; \endcode
    Matrix3D& operator=(const char* s)
    { 
      set_values(s); 
      return *this;
    }

    //! Sets all entries to zero (more efficient than assignement).
    void zero()
    {
      v_[0]=0.;  v_[1]=0.;  v_[2]=0.;
      v_[3]=0.;  v_[4]=0.;  v_[5]=0.;
      v_[6]=0.;  v_[7]=0.;  v_[8]=0.;
    }
    
    void identity()
    {
      v_[0]=1.;  v_[1]=0.;  v_[2]=0.;
      v_[3]=0.;  v_[4]=1.;  v_[5]=0.;
      v_[6]=0.;  v_[7]=0.;  v_[8]=1.;
    }
      
     
    //! Sets column j (0, 1 or 2) to Vector3D c.
    void set_column(int j, const Vector3D& c)
    {
      v_[0+j]=c[0];
      v_[3+j]=c[1];
      v_[6+j]=c[2];
    }
    
    //! returns the column length -- i is 0-based. 
    double column_length(int i) const 
    { return sqrt( v_[0+i]*v_[0+i] + v_[3+i]*v_[3+i] + v_[6+i]*v_[6+i] ); }


    double sub_det( int r, int c ) const
    {
      int r1 = 3 * ((r + 1) % 3);
      int r2 = 3 * ((r + 2) % 3);
      int c1 =     ((c + 1) % 3);
      int c2 =     ((c + 2) % 3);
      return v_[r1+c1] * v_[r2+c2] - v_[r2+c1] * v_[r1+c2];
    }
    
    // Matrix Operators
    friend bool operator==(const Matrix3D &lhs, const Matrix3D &rhs);
    friend bool operator!=(const Matrix3D &lhs, const Matrix3D &rhs);
    friend Matrix3D operator-( const Matrix3D& A );
    friend double Frobenius_2(const Matrix3D &A);
    friend Matrix3D transpose(const Matrix3D &A);
    inline Matrix3D& transpose();
    friend const Matrix3D operator+(const Matrix3D &A, const Matrix3D &B);
    friend const Matrix3D operator-(const Matrix3D &A, const Matrix3D &B) ;
    friend const Matrix3D operator*(const Matrix3D &A, const Matrix3D &B);
    inline Matrix3D& equal_mult_elem( const Matrix3D& A );
    friend const Matrix3D mult_element(const Matrix3D &A, const Matrix3D &B);
    inline Matrix3D& assign_product( const Matrix3D& A, const Matrix3D& B );
    friend void matmult(Matrix3D& C, const Matrix3D  &A, const Matrix3D &B);
    friend const Vector3D operator*(const Matrix3D  &A, const Vector3D &x);
    friend const Vector3D operator*(const Vector3D &x, const Matrix3D  &A);
    const Matrix3D operator*(double s) const;
    friend const Matrix3D operator*(double s, const Matrix3D &A);    
    void operator+=(const Matrix3D &rhs);
    void operator+=(const SymMatrix3D &rhs);
    void operator-=(const Matrix3D &rhs);
    void operator-=(const SymMatrix3D &rhs);
    void operator*=(double s);
    friend Matrix3D plus_transpose(const Matrix3D& A, const Matrix3D &B);
    Matrix3D& plus_transpose_equal(const Matrix3D &B);
    Matrix3D& outer_product(const Vector3D &v1, const Vector3D &v2);
    void fill_lower_triangle();

    //! \f$ v = A*x \f$
    friend void eqAx(Vector3D& v, const Matrix3D& A, const Vector3D& x);
    //! \f$ v += A*x \f$
    friend void plusEqAx(Vector3D& v, const Matrix3D& A, const Vector3D& x);
    friend void eqTransAx(Vector3D& v, const Matrix3D& A, const Vector3D& x); 
    //! \f$ v += A^T*x \f$
    friend void plusEqTransAx(Vector3D& v, const Matrix3D& A, const Vector3D& x);
     
    //! \f$ B += a*A \f$
    friend void plusEqaA(Matrix3D& B, const double a, const Matrix3D &A);

    //! determinant of matrix A, det(A).
    friend double det(const Matrix3D &A);

    //! \f$ B = A^{-1} \f$
    friend void inv(Matrix3D& B, const Matrix3D &A);

    //! \f$ B *= A^{-1} \f$
    friend void timesInvA(Matrix3D& B, const Matrix3D &A);

    //! \f$ Q*R = A \f$
    friend void QR(Matrix3D &Q, Matrix3D &R, const Matrix3D &A);

    size_t num_rows() const { return 3; }
    size_t num_cols() const { return 3; }

    //! returns a pointer to a row.
    inline double* operator[](unsigned i)
    {
      return v_ + 3*i;
    }

    //! returns a pointer to a row.
    inline const double* operator[](unsigned i) const
    {
      return v_ + 3*i;
    }
    
    inline double& operator()(unsigned short r, unsigned short c)
    {
      return v_[3*r+c];
    }
    inline double operator()(unsigned short r, unsigned short c) const
    {
      return v_[3*r+c];
    }
    
    
    
    inline Vector3D row(unsigned r) const
    {
      return Vector3D( v_ + 3*r );
    }
    
    inline Vector3D column( unsigned c ) const
    {
      return Vector3D( v_[c], v_[c+3], v_[c+6] );
    }

    inline bool positive_definite() const;
    
    inline SymMatrix3D upper() const
    {
      return SymMatrix3D( v_[0], v_[1], v_[2],
                                 v_[4], v_[5],
                                        v_[8] );
    }
    inline SymMatrix3D lower() const
    {
      return SymMatrix3D( v_[0], v_[3], v_[6],
                                 v_[4], v_[7],
                                        v_[8] );
    }
  };


  /* ***********  I/O  **************/

  inline std::ostream& operator<<(std::ostream &s, const Matrix3D &A)
  {
    for (size_t i=0; i<3; ++i)
      {
        for (size_t j=0; j<3; ++j)
          s << A[i][j] << " ";
        s << "\n";
      }
    return s;
  }

  inline std::istream& operator>>(std::istream &s, Matrix3D &A)
  {
    for (size_t i=0; i<3; i++)
      for (size_t j=0; j<3; j++)
        {
          s >>  A[i][j];
        }
    return s;
  }

  void Matrix3D::set_values(const char *s)
  {
    std::istringstream ins(s);
    ins >> *this;
  }

  // *********** matrix operators *******************

  // comparison functions
  inline bool operator==(const Matrix3D &lhs, const Matrix3D &rhs)
  {
    return lhs.v_[0] == rhs.v_[0]
        && lhs.v_[1] == rhs.v_[1]
        && lhs.v_[2] == rhs.v_[2]
        && lhs.v_[3] == rhs.v_[3]
        && lhs.v_[4] == rhs.v_[4]
        && lhs.v_[5] == rhs.v_[5]
        && lhs.v_[6] == rhs.v_[6]
        && lhs.v_[7] == rhs.v_[7]
        && lhs.v_[8] == rhs.v_[8];
  }
  inline bool operator!=(const Matrix3D &lhs, const Matrix3D &rhs)
  { return !(lhs == rhs); }

  inline Matrix3D operator-( const Matrix3D& A )
  {
    return Matrix3D( -A.v_[0],
                     -A.v_[1],
                     -A.v_[2],
                     -A.v_[3],
                     -A.v_[4],
                     -A.v_[5],
                     -A.v_[6],
                     -A.v_[7],
                     -A.v_[8] );
  }
                     
  //! \return A+B
  inline const Matrix3D operator+(const Matrix3D &A, 
                            const Matrix3D &B)
  {
    Matrix3D tmp(A);
    tmp += B;
    return tmp;
  }
  
  inline Matrix3D operator+( const Matrix3D& A, const SymMatrix3D& B )
  {
    return Matrix3D( A(0,0) + B[SymMatrix3D::T00],
                     A(0,1) + B[SymMatrix3D::T01],
                     A(0,2) + B[SymMatrix3D::T02],
                     A(1,0) + B[SymMatrix3D::T10],
                     A(1,1) + B[SymMatrix3D::T11],
                     A(1,2) + B[SymMatrix3D::T12],
                     A(2,0) + B[SymMatrix3D::T20],
                     A(2,1) + B[SymMatrix3D::T21],
                     A(2,2) + B[SymMatrix3D::T22] );
  }
  inline Matrix3D operator+( const SymMatrix3D& B, const Matrix3D& A )
    { return A + B; }

  //! \return A-B
  inline const Matrix3D operator-(const Matrix3D &A, 
                            const Matrix3D &B)
  {
    Matrix3D tmp(A);
    tmp -= B;
    return tmp;
  }
  
  inline Matrix3D operator-( const Matrix3D& A, const SymMatrix3D& B )
  {
    return Matrix3D( A(0,0) - B[SymMatrix3D::T00],
                     A(0,1) - B[SymMatrix3D::T01],
                     A(0,2) - B[SymMatrix3D::T02],
                     A(1,0) - B[SymMatrix3D::T10],
                     A(1,1) - B[SymMatrix3D::T11],
                     A(1,2) - B[SymMatrix3D::T12],
                     A(2,0) - B[SymMatrix3D::T20],
                     A(2,1) - B[SymMatrix3D::T21],
                     A(2,2) - B[SymMatrix3D::T22] );
  }
  inline Matrix3D operator-( const SymMatrix3D& B, const Matrix3D& A )
  {
    return Matrix3D( B[SymMatrix3D::T00] - A(0,0),
                     B[SymMatrix3D::T01] - A(0,1),
                     B[SymMatrix3D::T02] - A(0,2),
                     B[SymMatrix3D::T10] - A(1,0),
                     B[SymMatrix3D::T11] - A(1,1),
                     B[SymMatrix3D::T12] - A(1,2),
                     B[SymMatrix3D::T20] - A(2,0),
                     B[SymMatrix3D::T21] - A(2,1),
                     B[SymMatrix3D::T22] - A(2,2) );
  }
  
  inline Matrix3D& Matrix3D::equal_mult_elem( const Matrix3D& A )
  {
  	v_[0] *= A.v_[0];    
  	v_[1] *= A.v_[1];    
  	v_[2] *= A.v_[2];    
  	v_[3] *= A.v_[3];    
  	v_[4] *= A.v_[4];    
  	v_[5] *= A.v_[5];    
  	v_[6] *= A.v_[6];    
  	v_[7] *= A.v_[7];    
  	v_[8] *= A.v_[8];
    return *this;
  } 

    //! Multiplies entry by entry. This is NOT a matrix multiplication. 
  inline const Matrix3D mult_element(const Matrix3D &A, 
                               const Matrix3D &B)
  {
    Matrix3D tmp(A);
    tmp.equal_mult_elem(B);
    return tmp;
  }

  //! Return the square of the Frobenius norm of A, i.e. sum (diag (A' * A))
  inline double Frobenius_2(const Matrix3D &A)
  {
    return A.v_[0] * A.v_[0]
         + A.v_[1] * A.v_[1]
         + A.v_[2] * A.v_[2]
         + A.v_[3] * A.v_[3]
         + A.v_[4] * A.v_[4]
         + A.v_[5] * A.v_[5]
         + A.v_[6] * A.v_[6]
         + A.v_[7] * A.v_[7]
         + A.v_[8] * A.v_[8];
  }
  
  inline Matrix3D& Matrix3D::transpose()
  {
    double t;
    t = v_[1]; v_[1] = v_[3]; v_[3] = t;
    t = v_[2]; v_[2] = v_[6]; v_[6] = t;
    t = v_[5]; v_[5] = v_[7]; v_[7] = t;
    return *this;
  }
  
  inline Matrix3D transpose(const Matrix3D &A)
  {
    Matrix3D S;
//     size_t i;
//     for (i=0; i<3; ++i) {
//         S[size_t(0)][i] = A[i][0];
//         S[size_t(1)][i] = A[i][1];
//         S[size_t(2)][i] = A[i][2];
//     }
    S.v_[0]=A.v_[0]; S.v_[1]=A.v_[3]; S.v_[2]=A.v_[6];
    S.v_[3]=A.v_[1]; S.v_[4]=A.v_[4]; S.v_[5]=A.v_[7];
    S.v_[6]=A.v_[2]; S.v_[7]=A.v_[5]; S.v_[8]=A.v_[8];
    
    return S;
  }

  inline void Matrix3D::operator+=(const Matrix3D &rhs)
  {
      v_[0] += rhs.v_[0]; v_[1] += rhs.v_[1]; v_[2] += rhs.v_[2];
      v_[3] += rhs.v_[3]; v_[4] += rhs.v_[4]; v_[5] += rhs.v_[5];
      v_[6] += rhs.v_[6]; v_[7] += rhs.v_[7]; v_[8] += rhs.v_[8];
  }

  inline void Matrix3D::operator+=(const SymMatrix3D &rhs)
  {
      v_[0] += rhs[0]; v_[1] += rhs[1]; v_[2] += rhs[2];
      v_[3] += rhs[1]; v_[4] += rhs[3]; v_[5] += rhs[4];
      v_[6] += rhs[2]; v_[7] += rhs[4]; v_[8] += rhs[5];
  }

  inline void Matrix3D::operator-=(const Matrix3D &rhs)
  {
      v_[0] -= rhs.v_[0]; v_[1] -= rhs.v_[1]; v_[2] -= rhs.v_[2];
      v_[3] -= rhs.v_[3]; v_[4] -= rhs.v_[4]; v_[5] -= rhs.v_[5];
      v_[6] -= rhs.v_[6]; v_[7] -= rhs.v_[7]; v_[8] -= rhs.v_[8];
  }

  inline void Matrix3D::operator-=(const SymMatrix3D &rhs)
  {
      v_[0] -= rhs[0]; v_[1] -= rhs[1]; v_[2] -= rhs[2];
      v_[3] -= rhs[1]; v_[4] -= rhs[3]; v_[5] -= rhs[4];
      v_[6] -= rhs[2]; v_[7] -= rhs[4]; v_[8] -= rhs[5];
  }

  //! multiplies each entry by the scalar s
  inline void Matrix3D::operator*=(double s)
  {
      v_[0] *= s; v_[1] *= s; v_[2] *= s;
      v_[3] *= s; v_[4] *= s; v_[5] *= s;
      v_[6] *= s; v_[7] *= s; v_[8] *= s;
  }
  
  //! \f$ += B^T  \f$
  inline Matrix3D& Matrix3D::plus_transpose_equal( const Matrix3D& b )
  {
    if (&b == this) {
      v_[0] *= 2.0;
      v_[1] += v_[3];
      v_[2] += v_[6];
      v_[3]  = v_[1];
      v_[4] *= 2.0;
      v_[5] += v_[7];
      v_[6]  = v_[2];
      v_[7]  = v_[5];
      v_[8] *= 2.0;
    }
    else {
      v_[0] += b.v_[0];
      v_[1] += b.v_[3];
      v_[2] += b.v_[6];

      v_[3] += b.v_[1];
      v_[4] += b.v_[4];
      v_[5] += b.v_[7];

      v_[6] += b.v_[2];
      v_[7] += b.v_[5];
      v_[8] += b.v_[8];
    }
    return *this;
  }

  //! \f$ + B^T  \f$
  inline Matrix3D plus_transpose(const Matrix3D& A, const Matrix3D &B) 
  {
    Matrix3D tmp(A);
    tmp.plus_transpose_equal( B );
    return tmp;
  }

  //! Computes \f$ A = v_1 v_2^T \f$
  inline Matrix3D& Matrix3D::outer_product(const Vector3D  &v1, const Vector3D &v2)
  {
    // remember, matrix entries are v_[0] to v_[8].
    
    // diagonal
    v_[0] = v1[0]*v2[0];
    v_[4] = v1[1]*v2[1];
    v_[8] = v1[2]*v2[2];

    // upper triangular part
    v_[1] = v1[0]*v2[1];
    v_[2] = v1[0]*v2[2];
    v_[5] = v1[1]*v2[2];

    // lower triangular part
    v_[3] = v2[0]*v1[1];
    v_[6] = v2[0]*v1[2];
    v_[7] = v2[1]*v1[2];

    return *this;
  }

  inline void Matrix3D::fill_lower_triangle()
  {
    v_[3] = v_[1];
    v_[6] = v_[2];
    v_[7] = v_[5];
  } 

  //! \return A*B
  inline const Matrix3D operator*(const Matrix3D  &A, 
                            const Matrix3D &B)
  {
    Matrix3D tmp;
    tmp.assign_product( A, B );
    return tmp;
  }
  
  inline const Matrix3D operator*( const Matrix3D& A, 
                                   const SymMatrix3D& B )
  {
    return Matrix3D( A(0,0)*B[0] + A(0,1)*B[1] + A(0,2)*B[2],
                     A(0,0)*B[1] + A(0,1)*B[3] + A(0,2)*B[4],
                     A(0,0)*B[2] + A(0,1)*B[4] + A(0,2)*B[5],
                     
                     A(1,0)*B[0] + A(1,1)*B[1] + A(1,2)*B[2],
                     A(1,0)*B[1] + A(1,1)*B[3] + A(1,2)*B[4],
                     A(1,0)*B[2] + A(1,1)*B[4] + A(1,2)*B[5],
                     
                     A(2,0)*B[0] + A(2,1)*B[1] + A(2,2)*B[2],
                     A(2,0)*B[1] + A(2,1)*B[3] + A(2,2)*B[4],
                     A(2,0)*B[2] + A(2,1)*B[4] + A(2,2)*B[5] );
  }
  
  inline const Matrix3D operator*( const SymMatrix3D& B,
                                   const Matrix3D& A )
  {
    return Matrix3D( A(0,0)*B[0] + A(1,0)*B[1] + A(2,0)*B[2],
                     A(0,1)*B[0] + A(1,1)*B[1] + A(2,1)*B[2],
                     A(0,2)*B[0] + A(1,2)*B[1] + A(2,2)*B[2],

                     A(0,0)*B[1] + A(1,0)*B[3] + A(2,0)*B[4],
                     A(0,1)*B[1] + A(1,1)*B[3] + A(2,1)*B[4],
                     A(0,2)*B[1] + A(1,2)*B[3] + A(2,2)*B[4],

                     A(0,0)*B[2] + A(1,0)*B[4] + A(2,0)*B[5],
                     A(0,1)*B[2] + A(1,1)*B[4] + A(2,1)*B[5],
                     A(0,2)*B[2] + A(1,2)*B[4] + A(2,2)*B[5] );
  }
  
  inline const Matrix3D operator*( const SymMatrix3D& a,
                                   const SymMatrix3D& b )
  {
    return Matrix3D( a[0]*b[0] + a[1]*b[1] + a[2]*b[2],
                     a[0]*b[1] + a[1]*b[3] + a[2]*b[4],
                     a[0]*b[2] + a[1]*b[4] + a[2]*b[5],
                     
                     a[1]*b[0] + a[3]*b[1] + a[4]*b[2],
                     a[1]*b[1] + a[3]*b[3] + a[4]*b[4],
                     a[1]*b[2] + a[3]*b[4] + a[4]*b[5],
                     
                     a[2]*b[0] + a[4]*b[1] + a[5]*b[2],
                     a[2]*b[1] + a[4]*b[3] + a[5]*b[4],
                     a[2]*b[2] + a[4]*b[4] + a[5]*b[5] );
  }
                     
                     
   
   //! multiplies each entry by the scalar s
  inline const Matrix3D Matrix3D::operator*(double s) const
  {
    Matrix3D temp(*this);
    temp *= s;
    return temp;
  }
     //!friend function to allow for commutatative property of
     //! scalar mulitplication.
   inline const Matrix3D operator*(double s, const Matrix3D &A)
   {
     return (A.operator*(s));
   }
   

  inline Matrix3D& Matrix3D::assign_product( const Matrix3D& A, const Matrix3D& B)
  {
    v_[0] = A.v_[0]*B.v_[0] + A.v_[1]*B.v_[3] + A.v_[2]*B.v_[6];
    v_[1] = A.v_[0]*B.v_[1] + A.v_[1]*B.v_[4] + A.v_[2]*B.v_[7];
    v_[2] = A.v_[0]*B.v_[2] + A.v_[1]*B.v_[5] + A.v_[2]*B.v_[8];
    v_[3] = A.v_[3]*B.v_[0] + A.v_[4]*B.v_[3] + A.v_[5]*B.v_[6];
    v_[4] = A.v_[3]*B.v_[1] + A.v_[4]*B.v_[4] + A.v_[5]*B.v_[7];
    v_[5] = A.v_[3]*B.v_[2] + A.v_[4]*B.v_[5] + A.v_[5]*B.v_[8];
    v_[6] = A.v_[6]*B.v_[0] + A.v_[7]*B.v_[3] + A.v_[8]*B.v_[6];
    v_[7] = A.v_[6]*B.v_[1] + A.v_[7]*B.v_[4] + A.v_[8]*B.v_[7];
    v_[8] = A.v_[6]*B.v_[2] + A.v_[7]*B.v_[5] + A.v_[8]*B.v_[8];
    return *this;
  }
   

  //! \f$ C = A \times B \f$
  inline void matmult(Matrix3D& C, const Matrix3D  &A, const Matrix3D &B)
  {
    C.assign_product( A, B );
  }

  /*! \brief Computes \f$ A v \f$ . */
  inline const Vector3D operator*(const Matrix3D  &A, const Vector3D &x)
  {
    Vector3D tmp;
    eqAx( tmp, A, x );
    return tmp;
  }

  /*! \brief Computes \f$ v^T A \f$ .
    
      This function implicitly considers the transpose of vector x times
      the matrix A and it is implicit that the returned vector must be
      transposed. */
  inline const Vector3D operator*(const Vector3D &x, const Matrix3D  &A)
  {
    Vector3D tmp;
    eqTransAx( tmp, A, x );
    return tmp;
  }
   
  inline void eqAx(Vector3D& v, const Matrix3D& A, const Vector3D& x)
  {
     v.mCoords[0] = A.v_[0]*x[0] + A.v_[1]*x.mCoords[1] + A.v_[2]*x.mCoords[2];
     v.mCoords[1] = A.v_[3]*x[0] + A.v_[4]*x.mCoords[1] + A.v_[5]*x.mCoords[2];
     v.mCoords[2] = A.v_[6]*x[0] + A.v_[7]*x.mCoords[1] + A.v_[8]*x.mCoords[2];
  }
   
  inline void plusEqAx(Vector3D& v, const Matrix3D& A, const Vector3D& x)
  {
     v.mCoords[0] += A.v_[0]*x[0] + A.v_[1]*x.mCoords[1] + A.v_[2]*x.mCoords[2];
     v.mCoords[1] += A.v_[3]*x[0] + A.v_[4]*x.mCoords[1] + A.v_[5]*x.mCoords[2];
     v.mCoords[2] += A.v_[6]*x[0] + A.v_[7]*x.mCoords[1] + A.v_[8]*x.mCoords[2];
  }
  
  inline void eqTransAx(Vector3D& v, const Matrix3D& A, const Vector3D& x)
  {
     v.mCoords[0] = A.v_[0]*x.mCoords[0] + A.v_[3]*x.mCoords[1] + A.v_[6]*x.mCoords[2];
     v.mCoords[1] = A.v_[1]*x.mCoords[0] + A.v_[4]*x.mCoords[1] + A.v_[7]*x.mCoords[2];
     v.mCoords[2] = A.v_[2]*x.mCoords[0] + A.v_[5]*x.mCoords[1] + A.v_[8]*x.mCoords[2];
  } 
   
  inline void plusEqTransAx(Vector3D& v, const Matrix3D& A, const Vector3D& x)
  {
     v.mCoords[0] += A.v_[0]*x.mCoords[0] + A.v_[3]*x.mCoords[1] + A.v_[6]*x.mCoords[2];
     v.mCoords[1] += A.v_[1]*x.mCoords[0] + A.v_[4]*x.mCoords[1] + A.v_[7]*x.mCoords[2];
     v.mCoords[2] += A.v_[2]*x.mCoords[0] + A.v_[5]*x.mCoords[1] + A.v_[8]*x.mCoords[2];
  }
   
  inline void plusEqaA(Matrix3D& B, const double a, const Matrix3D &A) {
    B.v_[0] += a*A.v_[0]; B.v_[1] += a*A.v_[1]; B.v_[2] += a*A.v_[2]; 
    B.v_[3] += a*A.v_[3]; B.v_[4] += a*A.v_[4]; B.v_[5] += a*A.v_[5];
    B.v_[6] += a*A.v_[6]; B.v_[7] += a*A.v_[7]; B.v_[8] += a*A.v_[8];
  }

  inline double det(const Matrix3D &A) {
    return (  A.v_[0]*(A.v_[4]*A.v_[8]-A.v_[7]*A.v_[5])
            -A.v_[1]*(A.v_[3]*A.v_[8]-A.v_[6]*A.v_[5])
            +A.v_[2]*(A.v_[3]*A.v_[7]-A.v_[6]*A.v_[4]) );
  }

  inline void inv(Matrix3D &Ainv, const Matrix3D &A) {
    double inv_detA = 1.0 / (det(A));
      //First row of Ainv
    Ainv.v_[0] = inv_detA*( A.v_[4]*A.v_[8]-A.v_[5]*A.v_[7] );
    Ainv.v_[1] = inv_detA*( A.v_[2]*A.v_[7]-A.v_[8]*A.v_[1] );
    Ainv.v_[2] = inv_detA*( A.v_[1]*A.v_[5]-A.v_[4]*A.v_[2] );
      //Second row of Ainv
    Ainv.v_[3] = inv_detA*( A.v_[5]*A.v_[6]-A.v_[8]*A.v_[3] );
    Ainv.v_[4] = inv_detA*( A.v_[0]*A.v_[8]-A.v_[6]*A.v_[2] );
    Ainv.v_[5] = inv_detA*( A.v_[2]*A.v_[3]-A.v_[5]*A.v_[0] );
      //Third row of Ainv
    Ainv.v_[6] = inv_detA*( A.v_[3]*A.v_[7]-A.v_[6]*A.v_[4] );
    Ainv.v_[7] = inv_detA*( A.v_[1]*A.v_[6]-A.v_[7]*A.v_[0] );
    Ainv.v_[8] = inv_detA*( A.v_[0]*A.v_[4]-A.v_[3]*A.v_[1] );
  }

  inline void timesInvA(Matrix3D& B, const Matrix3D &A) {

    Matrix3D Ainv;
    inv( Ainv, A );
    B = B*Ainv;
  }

  inline void QR(Matrix3D &Q, Matrix3D &R, const Matrix3D &A) {
    // Compute the QR factorization of A.  This code uses the
    // Modified Gram-Schmidt method for computing the factorization.
    // The Householder version is more stable, but costs twice as many
    // floating point operations.

    Q = A;

    R[0][0] = sqrt(Q[0][0]*Q[0][0] + Q[1][0]*Q[1][0] + Q[2][0]*Q[2][0]);
    double temp_dbl = 1.0/R[0][0];
    R[1][0] = 0.0L;
    R[2][0] = 0.0L;
      //Q[0][0] /= R[0][0];
      //Q[1][0] /= R[0][0];
      //Q[2][0] /= R[0][0];
    Q[0][0] *= temp_dbl;
    Q[1][0] *= temp_dbl;
    Q[2][0] *= temp_dbl;
    

    R[0][1]  = Q[0][0]*Q[0][1] + Q[1][0]*Q[1][1] + Q[2][0]*Q[2][1];
    Q[0][1] -= Q[0][0]*R[0][1];
    Q[1][1] -= Q[1][0]*R[0][1];
    Q[2][1] -= Q[2][0]*R[0][1];

    R[0][2]  = Q[0][0]*Q[0][2] + Q[1][0]*Q[1][2] + Q[2][0]*Q[2][2];
    Q[0][2] -= Q[0][0]*R[0][2];
    Q[1][2] -= Q[1][0]*R[0][2];
    Q[2][2] -= Q[2][0]*R[0][2];

    R[1][1] = sqrt(Q[0][1]*Q[0][1] + Q[1][1]*Q[1][1] + Q[2][1]*Q[2][1]);
    temp_dbl = 1.0 / R[1][1];
    R[2][1] = 0.0L;
//     Q[0][1] /= R[1][1];
//     Q[1][1] /= R[1][1];
//     Q[2][1] /= R[1][1];
    Q[0][1] *= temp_dbl;
    Q[1][1] *= temp_dbl;
    Q[2][1] *= temp_dbl;

    
    R[1][2]  = Q[0][1]*Q[0][2] + Q[1][1]*Q[1][2] + Q[2][1]*Q[2][2];
    Q[0][2] -= Q[0][1]*R[1][2];
    Q[1][2] -= Q[1][1]*R[1][2];
    Q[2][2] -= Q[2][1]*R[1][2];
  
    R[2][2] = sqrt(Q[0][2]*Q[0][2] + Q[1][2]*Q[1][2] + Q[2][2]*Q[2][2]);
    temp_dbl = 1.0 / R[2][2];
    
//     Q[0][2] /= R[2][2];
//     Q[1][2] /= R[2][2];
//     Q[2][2] /= R[2][2];
    Q[0][2] *= temp_dbl;
    Q[1][2] *= temp_dbl;
    Q[2][2] *= temp_dbl;
    
    return;
  }
  
inline bool Matrix3D::positive_definite() const
{
  // A = B + C
 //where
 // B = (A + transpose(A))/2
 // C = (A - transpose(A))/2
 // B is always a symmetric matrix and 
 // A is positive definite iff B is positive definite.
 Matrix3D B(*this);
 B.plus_transpose_equal( *this );
 B *= 0.5;
 
 // Sylvester's Criterion for positive definite symmetric matrix
 return (B[0][0] > 0.0) && (B.sub_det(2,2) > 0.0) && (det(B) > 0.0);
}

} // namespace Mesquite

#endif // Matrix3D_hpp
