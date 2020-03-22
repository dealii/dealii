/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file MsqMatrix.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_MSQ_MATRIX_HPP
#define MSQ_MSQ_MATRIX_HPP

#include "Mesquite.hpp"

#include <string>
#include <istream>
#include <ostream>
#include <sstream>

namespace MESQUITE_NS {

/**\brief Fixed-size matrix class
 *
 * This class implements a fixed-size 2-dimensional matrix.
 * The actual size is specified with template parameters.
 */
template <unsigned R, unsigned C> 
class MsqMatrix
{
protected:
  double m[R*C];

public:
  typedef MsqMatrix<R,C> my_type;
  
  enum { ROWS = R, COLS = C };

    /** Constructor for uninitialized matrix */
  MsqMatrix()                                         { }
    /** Initialize diagonal values, zero others */
  MsqMatrix( double v )                               { diag(v); }
    /** Initialize to an array of values */
  MsqMatrix( const double* v )                        { set(v); }
    /** Initialize from 2D array */
  MsqMatrix( const double v[R][C] )                   { set(v); }
    /** Initialize with column vectors */
  MsqMatrix( const MsqMatrix<R,1>* c )                { set_columns(c); }
    /** Initialize with row vectors */
  MsqMatrix( const MsqMatrix<1,C>* r )                { set_rows(r); }
    /** Parse initial values from string */
  MsqMatrix( const char* s )                          { set(s); }
    /** Parse initial values from string */
  MsqMatrix( const std::string& s )               { set(s); }
    /** Initialize to the minor of a larger matrix
     *  This matrix is the passed matrix with the 
     *  specified row and column removed.
     */
  MsqMatrix( const MsqMatrix<R+1,C+1>& m, unsigned r, unsigned c )
                                                      { make_minor(m,r,c); }
  
  MsqMatrix<R,C>& operator=( double v )                 { set(v); return *this; }
  MsqMatrix<R,C>& operator=( const double* v )          { set(v); return *this; }
  MsqMatrix<R,C>& operator=( const char* s )            { set(s); return *this; }
  MsqMatrix<R,C>& operator=( const std::string& s ) { set(s); return *this; }
  
  double& operator()( unsigned r, unsigned c )        { return m[r*C+c]; }
  double  operator()( unsigned r, unsigned c ) const  { return m[r*C+c]; }
  double* data()                                      { return m; }
  const double* data() const                          { return m; }
  
  void zero()                                         { set(0.0);  }
  void identity()                                     { diag(1.0); }
  void set( double v )        { for (unsigned i = 0; i < R*C; ++i) m[i] = v; }
  void set( const double* v ) { for (unsigned i = 0; i < R*C; ++i) m[i] = v[i]; }
  void set( const double v[R][C] );
  void set( const char* s )   { std::istringstream i(s); i >> *this; }
  void set( const std::string& s ) { set( s.c_str() ); }
    /** Set diagonal value to passed values, others to zero. */
  inline void diag( double v );
    /** Set diagonal values to passed values, others to zero. */
  inline void diag( const double* v );
    /** Set this matrix to the minor of a larger matrix */
  inline void make_minor( const MsqMatrix<R+1,C+1>& m, unsigned r, unsigned c );
  
    /** *this += transpose(other) */
  inline MsqMatrix<R,C>& assign_add_transpose( const MsqMatrix<C,R>& other );
    /** *this = s*m */
  inline MsqMatrix<R,C>& assign_product( double s, const MsqMatrix<R,C>& m );
    /** *this += s*m */
  inline MsqMatrix<R,C>& assign_add_product( double s, const MsqMatrix<R,C>& m );
    /** multiply each element by the cooresponding element in m */
  inline MsqMatrix<R,C>& assign_multiply_elements( const MsqMatrix<R,C>& m );
  
  inline MsqMatrix<R,C>& operator+=( const MsqMatrix<R,C>& other );
  inline MsqMatrix<R,C>& operator-=( const MsqMatrix<R,C>& other );
  inline MsqMatrix<R,C>& operator+=( double scalar );
  inline MsqMatrix<R,C>& operator-=( double scalar );
  inline MsqMatrix<R,C>& operator*=( double scalar );
  inline MsqMatrix<R,C>& operator/=( double scalar );
  
  //      MsqMatrix<1,C>& row( unsigned r )       { return *(MsqMatrix<1,C>*)(m+C*r); }
  //const MsqMatrix<1,C>& row( unsigned r ) const { return *(MsqMatrix<1,C>*)(m+C*r); }
  MsqMatrix<1,C> row( unsigned r ) const { return MsqMatrix<1,C>(m+C*r); }
  MsqMatrix<R,1> column( unsigned c ) const;
  MsqMatrix<1,R> column_transpose( unsigned c ) const;
  
  void set_row( unsigned r, const MsqMatrix<1,C>& v );
  void add_row( unsigned r, const MsqMatrix<1,C>& v );
  void set_row_transpose( unsigned r, const MsqMatrix<C,1>& v );
  void set_rows( const MsqMatrix<1,C>* v );
  void set_column( unsigned c, const MsqMatrix<R,1>& v );
  void add_column( unsigned c, const MsqMatrix<R,1>& v );
  void set_column_transpose( unsigned c, const MsqMatrix<1,R>& v );
  void set_columns( const MsqMatrix<R,1>* v );
};

template <> 
class MsqMatrix<1,1>
{
protected:
  double m;

public:
  typedef MsqMatrix<1,1> my_type;
  
  enum { ROWS = 1, COLS = 1 };

    /** Constructor for uninitialized matrix */
  MsqMatrix()                                         { }
    /** Initialize diagonal values, zero others */
  MsqMatrix( double v )                               : m(v) {}
    /** Initialize to an array of values */
  MsqMatrix( const double* v )                        : m(*v) {}
    /** Parse initial values from string */
  MsqMatrix( const char* s )                          { set(s); }
    /** Parse initial values from string */
  MsqMatrix( const std::string& s )               { set(s); }
    /** Initialize to the minor of a larger matrix
     *  This matrix is the passed matrix with the 
     *  specified row and column removed.
     */
  MsqMatrix( const MsqMatrix<2,2>& M, unsigned r, unsigned c ) :  m(M(r,c)) {}
  
  MsqMatrix<1,1>& operator=( double v )                 { m = v; return *this; }
  MsqMatrix<1,1>& operator=( const double* v )          { m = *v; return *this; }
  MsqMatrix<1,1>& operator=( const char* s )            { set(s); return *this; }
  MsqMatrix<1,1>& operator=( const std::string& s ) { set(s); return *this; }
  
  double& operator()( unsigned, unsigned )        { return m; }
  double  operator()( unsigned, unsigned ) const  { return m; }
  double* data()                                      { return &m; }
  const double* data() const                          { return &m; }
  
  void zero()                                         { m = 0.0; }
  void identity()                                     { m = 1.0; }
  void set( double v )        { m = v; }
  void set( const double* v ) { m= *v; }
  void set( const char* s )   { std::istringstream i(s); i >> m; }
  void set( const std::string& s ) { set( s.c_str() ); }
    /** Set diagonal value to passed values, others to zero. */
  inline void diag( double v ) { m = v; }
    /** Set diagonal values to passed values, others to zero. */
  inline void diag( const double* v ) { m = *v; }
    /** Set this matrix to the minor of a larger matrix */
  inline void make_minor( const MsqMatrix<2,2>& M, unsigned r, unsigned c )
    { m = M(r,c); }
  
    /** *this += transpose(other) */
  inline MsqMatrix<1,1>& assign_add_transpose( const MsqMatrix<1,1>& other )
    { m += other.m; return *this; }
    /** *this = s*m */
  inline MsqMatrix<1,1>& assign_product( double s, const MsqMatrix<1,1>& other )
    { m = s*other.m; return *this; }
    /** *this += s*m */
  inline MsqMatrix<1,1>& assign_add_product( double s, const MsqMatrix<1,1>& other )
    { m += s*other.m; return *this; }
    /** multiply each element by the cooresponding element in m */
  inline MsqMatrix<1,1>& assign_multiply_elements( const MsqMatrix<1,1>& other )
    { m *= other.m; return *this; }
    
  operator double () const
    { return m; }
};

/** \brief Vector is a 1xL Matrix 
 *
 * Define a Vector as a 1xL Matrix
 * Add single-index access operators
 */
template <unsigned L> 
class MsqVector : public MsqMatrix<L,1>
{
public:
  MsqVector()                                         { }
  MsqVector( double v )                               { MsqMatrix<L,1>::set(v); }
  MsqVector( const double* v )                        { MsqMatrix<L,1>::set(v); }
  MsqVector( const char* s )                          { MsqMatrix<L,1>::set(s); }
  MsqVector( const std::string& s )               { MsqMatrix<L,1>::set(s); }
  MsqVector( const MsqMatrix<L,1>& m) : MsqMatrix<L,1>(m) {}
  
  double& operator[](unsigned idx)                    { return MsqMatrix<L,1>::operator()(idx,0); }
  double  operator[](unsigned idx) const              { return MsqMatrix<L,1>::operator()(idx,0); }
  double& operator()(unsigned idx)                    { return MsqMatrix<L,1>::operator()(idx,0); }
  double  operator()(unsigned idx) const              { return MsqMatrix<L,1>::operator()(idx,0); }

  MsqVector<L>& operator=( const MsqMatrix<L,1>& m )
    { MsqMatrix<L,1>::operator=(m); return *this; }
};

/**\brief A subclass of MsqMatrix that behaves more like the old Matrix3D class
 *
 * A MsqMartix that behaves more like the old Matrx3D class.
 * * Constructors are different--be careful.
 * * Adds operator[] for c-array style access.
 */
template <unsigned R, unsigned C>
class MsqMatrixA : public MsqMatrix<R,C>
{
public:
    /** Initialize all elements to zero */
  MsqMatrixA()                                         { MsqMatrix<R,C>::set(0); }
    /** Initialize all elements to passed value */
  MsqMatrixA( double v )                               { MsqMatrix<R,C>::set(v); }
  MsqMatrixA( const double* v )                        { MsqMatrix<R,C>::set(v); }
  MsqMatrixA( const char* s )                          { MsqMatrix<R,C>::set(s); }
  MsqMatrixA( const std::string& s )               { MsqMatrix<R,C>::set(s); }
  MsqMatrixA( const MsqMatrix<R,C>& m) : MsqMatrix<R,C>(m) {}
  MsqMatrixA( const MsqMatrix<R+1,C+1>& m, unsigned r, unsigned c )
                                                       { MsqMatrix<R,C>::make_minor(m,r,c); }
  
  double* operator[]( unsigned idx )                  { return MsqMatrix<R,C>::m+C*idx;  }
  const double* operator[]( unsigned idx ) const      { return MsqMatrix<R,C>::m+C*idx;  }

  MsqMatrixA<R,C>& operator=( const MsqMatrixA<R,C>& m )
    { MsqMatrixA<R,C>::operator=(m); return *this; }
};

template <unsigned R, unsigned C> inline
void MsqMatrix<R,C>::set( const double v[R][C] )
{
  for (unsigned r = 0; r < R; ++r)
    for (unsigned c = 0; c < C; ++c)
      operator()(r,c) = v[r][c];
}

template <unsigned R, unsigned C> inline
void MsqMatrix<R,C>::diag( double v ) {
  //for (unsigned r = 0; r < R; ++r)
  //  for (unsigned c = 0; c < C; ++c)
  //    operator()(r,c) = (r == c) ? v : 0.0;

  switch (R) {
    default: for (unsigned r = 4; r < R; ++r)
      switch (C) {
        default: for (unsigned k = 4; k < C; ++k)
                   operator()(r,k) = r == k ? v : 0.0;
        case 4:    operator()(r,3) = 0.0;
        case 3:    operator()(r,2) = 0.0;
        case 2:    operator()(r,1) = 0.0;
        case 1:    operator()(r,0) = 0.0;
      }
    case 4:
      switch (C) {
        default: for (unsigned k = 4; k < C; ++k)
                   operator()(3,k) = 0.0;
        case 4:    operator()(3,3) = v;
        case 3:    operator()(3,2) = 0.0;
        case 2:    operator()(3,1) = 0.0;
        case 1:    operator()(3,0) = 0.0;
      }
    case 3:
      switch (C) {
        default: for (unsigned k = 4; k < C; ++k)
                   operator()(2,k) = 0.0;
        case 4:    operator()(2,3) = 0.0;
        case 3:    operator()(2,2) = v;
        case 2:    operator()(2,1) = 0.0;
        case 1:    operator()(2,0) = 0.0;
      }
    case 2:
      switch (C) {
        default: for (unsigned k = 4; k < C; ++k)
                   operator()(1,k) = 0.0;
        case 4:    operator()(1,3) = 0.0;
        case 3:    operator()(1,2) = 0.0;
        case 2:    operator()(1,1) = v;
        case 1:    operator()(1,0) = 0.0;
      }
    case 1:
      switch (C) {
        default: for (unsigned k = 4; k < C; ++k)
                   operator()(0,k) = 0.0;
        case 4:    operator()(0,3) = 0.0;
        case 3:    operator()(0,2) = 0.0;
        case 2:    operator()(0,1) = 0.0;
        case 1:    operator()(0,0) = v;
      }
  }
}

template <unsigned R, unsigned C> inline
void MsqMatrix<R,C>::diag( const double* v ) {
  //for (unsigned r = 0; r < R; ++r)
  //  for (unsigned c = 0; c < C; ++c)
  //    operator()(r,c) = (r == c) ? v[r] : 0.0;

  switch (R) {
    default: for (unsigned r = 4; r < R; ++r)
      switch (C) {
        default: for (unsigned k = 4; k < C; ++k)
                   operator()(r,k) = r == k ? v[r] : 0.0;
        case 4:    operator()(r,3) = 0.0;
        case 3:    operator()(r,2) = 0.0;
        case 2:    operator()(r,1) = 0.0;
        case 1:    operator()(r,0) = 0.0;
      }
    case 4:
      switch (C) {
        default: for (unsigned k = 4; k < C; ++k)
                   operator()(3,k) = 0.0;
        case 4:    operator()(3,3) = v[3];
        case 3:    operator()(3,2) = 0.0;
        case 2:    operator()(3,1) = 0.0;
        case 1:    operator()(3,0) = 0.0;
      }
    case 3:
      switch (C) {
        default: for (unsigned k = 4; k < C; ++k)
                   operator()(2,k) = 0.0;
        case 4:    operator()(2,3) = 0.0;
        case 3:    operator()(2,2) = v[2];
        case 2:    operator()(2,1) = 0.0;
        case 1:    operator()(2,0) = 0.0;
      }
    case 2:
      switch (C) {
        default: for (unsigned k = 4; k < C; ++k)
                   operator()(1,k) = 0.0;
        case 4:    operator()(1,3) = 0.0;
        case 3:    operator()(1,2) = 0.0;
        case 2:    operator()(1,1) = v[1];
        case 1:    operator()(1,0) = 0.0;
      }
    case 1:
      switch (C) {
        default: for (unsigned k = 4; k < C; ++k)
                   operator()(0,k) = 0.0;
        case 4:    operator()(0,3) = 0.0;
        case 3:    operator()(0,2) = 0.0;
        case 2:    operator()(0,1) = 0.0;
        case 1:    operator()(0,0) = v[0];
      }
  }
}

/**\brief Extract minor of a matrix and assign to *this
 *
 * Given a matrix m, a row r and an column c, set *this to 
 * the matrix that is m with row r and column c deleted.
 */
template <unsigned R, unsigned C> inline
void MsqMatrix<R,C>::make_minor( const MsqMatrix<R+1,C+1>& M, unsigned r, unsigned c )
{
  for (unsigned i = 0; i < r; ++i) {
    for (unsigned j = 0; j < c; ++j)
      operator()(i,j) = M(i,j);
    for (unsigned j = c; j < C; ++j)
      operator()(i,j) = M(i,j+1);
  }
  for (unsigned i = r; i < R; ++i) {
    for (unsigned j = 0; j < c; ++j)
      operator()(i,j) = M(i+1,j);
    for (unsigned j = c; j < C; ++j)
      operator()(i,j) = M(i+1,j+1);
  }
}

template <unsigned R, unsigned C> inline
void MsqMatrix<R,C>::set_row( unsigned r, const MsqMatrix<1,C>& v )
  { for (unsigned i = 0; i < C; ++i) operator()(r,i) = v(0,i); }
template <unsigned R, unsigned C> inline
void MsqMatrix<R,C>::add_row( unsigned r, const MsqMatrix<1,C>& v )
  { for (unsigned i = 0; i < C; ++i) operator()(r,i) += v(0,i); }
template <unsigned R, unsigned C> inline
void MsqMatrix<R,C>::set_row_transpose( unsigned r, const MsqMatrix<C,1>& v )
  { for (unsigned i = 0; i < C; ++i) operator()(r,i) = v(i,0); }
template <unsigned R, unsigned C> inline
void MsqMatrix<R,C>::set_rows( const MsqMatrix<1,C>* v )
  { for( unsigned r = 0; r < R; ++r) set_row( r, v[r] ); }

template <unsigned R, unsigned C> inline
void MsqMatrix<R,C>::set_column( unsigned c, const MsqMatrix<R,1>& v )
  { for (unsigned i = 0; i < R; ++i) operator()(i,c) = v(i,0); }
template <unsigned R, unsigned C> inline
void MsqMatrix<R,C>::add_column( unsigned c, const MsqMatrix<R,1>& v )
  { for (unsigned i = 0; i < R; ++i) operator()(i,c) += v(i,0); }
template <unsigned R, unsigned C> inline
void MsqMatrix<R,C>::set_column_transpose( unsigned c, const MsqMatrix<1,R>& v )
  { for (unsigned i = 0; i < R; ++i) operator()(i,c) = v(0,i); }
template <unsigned R, unsigned C> inline
void MsqMatrix<R,C>::set_columns( const MsqMatrix<R,1>* v )
  { for( unsigned c = 0; c < C; ++c) set_column( c, v[c] ); }
  

template <unsigned R, unsigned C> inline
MsqMatrix<R,1> MsqMatrix<R,C>::column( unsigned c ) const
{
  MsqMatrix<R,1> result;
  for (unsigned i = 0; i < R; ++i)
    result(i,0) = operator()(i,c);
  return result;
}

template <unsigned R, unsigned C> inline
MsqMatrix<1,R> MsqMatrix<R,C>::column_transpose( unsigned c ) const
{
  MsqMatrix<1,R> result;
  for (unsigned i = 0; i < R; ++i)
    result(0,i) = operator()(i,c);
  return result;
}

/**\brief Set a subset of this matrix to some other matrix */
template <unsigned R1, unsigned C1, unsigned R2, unsigned C2> inline
void set_region( MsqMatrix<R1,C1>& d,
                 unsigned r, unsigned c, 
                 MsqMatrix<R2,C2>& s )
{
  const unsigned rmax = r+R2 > R1 ? R1 : r+R2;
  const unsigned cmax = c+C2 > C1 ? C1 : c+C2;
  for (unsigned i = r; i < rmax; ++i)
    for (unsigned j = c; j < cmax; ++j)
      d(i,j) = s(i-r,j-c);
}


template <unsigned R, unsigned C> inline
MsqMatrix<R,C>& MsqMatrix<R,C>::assign_add_transpose( const MsqMatrix<C,R>& other )
{
  for (unsigned r = 0; r < R; ++r)
    for (unsigned c = 0; c < C; ++c)
      operator()(r,c) += other(c,r);
  return *this;
}  

template <unsigned R, unsigned C> inline
MsqMatrix<R,C>& MsqMatrix<R,C>::assign_multiply_elements( const MsqMatrix<R,C>& other )
  { for (unsigned i = 0; i < R*C; ++i)  m[i] *= other.data()[i]; return *this; }  

template <unsigned R, unsigned C> inline
MsqMatrix<R,C>& MsqMatrix<R,C>::assign_product( double s, const MsqMatrix<R,C>& other )
  { for (unsigned i = 0; i < R*C; ++i)  m[i] = s * other.data()[i]; return *this; }  

template <unsigned R, unsigned C> inline
MsqMatrix<R,C>& MsqMatrix<R,C>::assign_add_product( double s, const MsqMatrix<R,C>& other )
  { for (unsigned i = 0; i < R*C; ++i) m[i] += s * other.data()[i]; return *this; }  

template <unsigned R, unsigned C> inline
MsqMatrix<R,C>& MsqMatrix<R,C>::operator+=( const MsqMatrix<R,C>& other )
  { for (unsigned i = 0; i < R*C; ++i)  m[i] += other.data()[i]; return *this; }  
  
template <unsigned R, unsigned C> inline
MsqMatrix<R,C>& MsqMatrix<R,C>::operator-=( const MsqMatrix<R,C>& other )
  { for (unsigned i = 0; i < R*C; ++i)  m[i] -= other.data()[i]; return *this; }  
  
template <unsigned R, unsigned C> inline
MsqMatrix<R,C>& MsqMatrix<R,C>::operator+=( double s )
  { for (unsigned i = 0; i < R*C; ++i)  m[i] += s; return *this; }  
  
template <unsigned R, unsigned C> inline
MsqMatrix<R,C>& MsqMatrix<R,C>::operator-=( double s )
  { for (unsigned i = 0; i < R*C; ++i)  m[i] -= s; return *this; }  
  
template <unsigned R, unsigned C> inline
MsqMatrix<R,C>& MsqMatrix<R,C>::operator*=( double s )
  { for (unsigned i = 0; i < R*C; ++i)  m[i] *= s; return *this; }  
  
template <unsigned R, unsigned C> inline
MsqMatrix<R,C>& MsqMatrix<R,C>::operator/=( double s )
  { for (unsigned i = 0; i < R*C; ++i)  m[i] /= s; return *this; }  


template <unsigned R, unsigned C> inline
MsqMatrix<R,C> operator-( const MsqMatrix<R,C>& m )
{
  MsqMatrix<R,C> result;
  for (unsigned i = 0; i < R*C; ++i)
    result.data()[i] = -m.data()[i];
  return result;
}

template <unsigned R, unsigned C> inline
MsqMatrix<R,C> operator+( const MsqMatrix<R,C>& m, double s )
  { MsqMatrix<R,C> tmp(m); tmp += s; return tmp; }

template <unsigned R, unsigned C> inline
MsqMatrix<R,C> operator+( double s, const MsqMatrix<R,C>& m )
  { MsqMatrix<R,C> tmp(m); tmp += s; return tmp; }

template <unsigned R, unsigned C> inline
MsqMatrix<R,C> operator+( const MsqMatrix<R,C>& A, const MsqMatrix<R,C>& B )
  { MsqMatrix<R,C> tmp(A); tmp += B; return tmp; }


template <unsigned R, unsigned C> inline
MsqMatrix<R,C> operator-( const MsqMatrix<R,C>& m, double s )
  { MsqMatrix<R,C> tmp(m); tmp -= s; return tmp; }

template <unsigned R, unsigned C> inline
MsqMatrix<R,C> operator-( double s, const MsqMatrix<R,C>& m )
  { MsqMatrix<R,C> tmp(m); tmp -= s; return tmp; }

template <unsigned R, unsigned C> inline
MsqMatrix<R,C> operator-( const MsqMatrix<R,C>& A, const MsqMatrix<R,C>& B )
  { MsqMatrix<R,C> tmp(A); tmp -= B; return tmp; }


template <unsigned R, unsigned C> inline
MsqMatrix<R,C> operator*( const MsqMatrix<R,C>& m, double s )
  { MsqMatrix<R,C> tmp(m); tmp *= s; return tmp; }

template <unsigned R, unsigned C> inline
MsqMatrix<R,C> operator*( double s, const MsqMatrix<R,C>& m )
  { MsqMatrix<R,C> tmp(m); tmp *= s; return tmp; }

template <unsigned R, unsigned RC, unsigned C> inline
double multiply_helper_result_val( unsigned r, unsigned c,
                                   const MsqMatrix<R,RC>& A,
                                   const MsqMatrix<RC,C>& B )
{
  double tmp = A(r,0)*B(0,c);
  switch (RC) {
    default: for (unsigned k = 6; k < RC; ++k)
            tmp += A(r,k)*B(k,c);
    case 6: tmp += A(r,5)*B(5,c);
    case 5: tmp += A(r,4)*B(4,c);
    case 4: tmp += A(r,3)*B(3,c);
    case 3: tmp += A(r,2)*B(2,c);
    case 2: tmp += A(r,1)*B(1,c);
    case 1: ;
  }
  return tmp;
}

template <unsigned R, unsigned RC, unsigned C> inline
MsqMatrix<R,C> operator*( const MsqMatrix<R,RC>& A, const MsqMatrix<RC,C>& B )
{
  //MsqMatrix<R,C> result(0.0);
  //for (unsigned i = 0; i < R; ++i)
  //  for (unsigned j = 0; j < C; ++j)
  //    for (unsigned k = 0; k < RC; ++k)
  //      result(i,j) += A(i,k) * B(k,j);
 
  MsqMatrix<R,C> result;
  for (unsigned r = 0; r < R; ++r)
    for (unsigned c = 0; c < C; ++c)
      result(r,c) = multiply_helper_result_val( r, c, A, B );
  
  return result;
}


template <unsigned R, unsigned C> inline
MsqMatrix<R,C> operator/( const MsqMatrix<R,C>& m, double s )
  { MsqMatrix<R,C> tmp(m); tmp /= s; return tmp; }

template <unsigned RC> inline
double cofactor( const MsqMatrix<RC,RC>& m, unsigned r, unsigned c )
{
  const double sign[] = { 1.0, -1.0 };
  return sign[(r+c)%2] * det( MsqMatrix<RC-1,RC-1>( m, r, c ) );
}

template <unsigned RC> inline
double det( const MsqMatrix<RC,RC>& m )
{
  double result = 0.0;
  for (unsigned i = 0; i < RC; ++i) 
    result += m(0,i) * cofactor<RC>( m, 0, i );
  return result;
}

//inline 
//double det( const MsqMatrix<1,1>& m )
//  { return m(0,0); }

inline 
double det( const MsqMatrix<2,2>& m )
  { return m(0,0)*m(1,1) - m(0,1)*m(1,0); }

inline
double det( const MsqMatrix<3,3>& m )
  { return m(0,0)*(m(1,1)*m(2,2) - m(2,1)*m(1,2))
         + m(0,1)*(m(2,0)*m(1,2) - m(1,0)*m(2,2))
         + m(0,2)*(m(1,0)*m(2,1) - m(2,0)*m(1,1)); }

inline
MsqMatrix<2,2> adj( const MsqMatrix<2,2>& m )
{
  MsqMatrix<2,2> result;
  result(0,0) =  m(1,1);
  result(0,1) = -m(0,1);
  result(1,0) = -m(1,0);
  result(1,1) =  m(0,0);
  return result;
}

inline
MsqMatrix<2,2> transpose_adj( const MsqMatrix<2,2>& m )
{
  MsqMatrix<2,2> result;
  result(0,0) =  m(1,1);
  result(1,0) = -m(0,1);
  result(0,1) = -m(1,0);
  result(1,1) =  m(0,0);
  return result;
}

inline
MsqMatrix<3,3> adj( const MsqMatrix<3,3>& m )
{
  MsqMatrix<3,3> result;
  
  result(0,0) = m(1,1)*m(2,2) - m(1,2)*m(2,1);
  result(0,1) = m(0,2)*m(2,1) - m(0,1)*m(2,2);
  result(0,2) = m(0,1)*m(1,2) - m(0,2)*m(1,1);
  
  result(1,0) = m(1,2)*m(2,0) - m(1,0)*m(2,2);
  result(1,1) = m(0,0)*m(2,2) - m(0,2)*m(2,0);
  result(1,2) = m(0,2)*m(1,0) - m(0,0)*m(1,2);
  
  result(2,0) = m(1,0)*m(2,1) - m(1,1)*m(2,0);
  result(2,1) = m(0,1)*m(2,0) - m(0,0)*m(2,1);
  result(2,2) = m(0,0)*m(1,1) - m(0,1)*m(1,0);
  
  return result;
}

inline
MsqMatrix<3,3> transpose_adj( const MsqMatrix<3,3>& m )
{
  MsqMatrix<3,3> result;
  
  result(0,0) = m(1,1)*m(2,2) - m(1,2)*m(2,1);
  result(0,1) = m(1,2)*m(2,0) - m(1,0)*m(2,2);
  result(0,2) = m(1,0)*m(2,1) - m(1,1)*m(2,0);
  
  result(1,0) = m(0,2)*m(2,1) - m(0,1)*m(2,2);
  result(1,1) = m(0,0)*m(2,2) - m(0,2)*m(2,0);
  result(1,2) = m(0,1)*m(2,0) - m(0,0)*m(2,1);
  
  result(2,0) = m(0,1)*m(1,2) - m(0,2)*m(1,1);
  result(2,1) = m(0,2)*m(1,0) - m(0,0)*m(1,2);
  result(2,2) = m(0,0)*m(1,1) - m(0,1)*m(1,0);
  
  return result;
}

template <unsigned R, unsigned C> inline
MsqMatrix<R,C> inverse( const MsqMatrix<R,C>& m )
{
  return adj(m) * (1.0 / det(m));
}

template <unsigned R, unsigned C> inline
MsqMatrix<R,C> transpose( const MsqMatrix<C,R>& m )
{
  MsqMatrix<R,C> result;
  for (unsigned r = 0; r < R; ++r)
    for (unsigned c = 0; c < C; ++c)
      result(r,c) = m(c,r);
  return result;
}
/*
template <unsigned R> inline
const MsqMatrix<R,1>& transpose( const MsqMatrix<1,R>& m )
  { return *reinterpret_cast<const MsqMatrix<R,1>*>(&m); }

template <unsigned C> inline
const MsqMatrix<1,C>& transpose( const MsqMatrix<C,1>& m )
  { return *reinterpret_cast<const MsqMatrix<1,C>*>(&m); }
*/
template <unsigned RC> inline
double trace( const MsqMatrix<RC,RC>& m )
{
  double result = m(0,0);
  for (unsigned i = 1; i < RC; ++i)
    result += m(i,i);
  return result;
}

template <unsigned R, unsigned C> inline
double sqr_Frobenius( const MsqMatrix<R,C>& m )
{
  double result = *m.data() * *m.data();
  for (unsigned i = 1; i < R*C; ++i)
    result += m.data()[i] * m.data()[i];
  return result;
}

template <unsigned R, unsigned C> inline
double Frobenius( const MsqMatrix<R,C>& m )
  { return std::sqrt( sqr_Frobenius<R,C>( m ) ); }

template <unsigned R, unsigned C> inline
bool operator==( const MsqMatrix<R,C>& A, const MsqMatrix<R,C>& B ) 
{
  for (unsigned i = 0; i < R*C; ++i)
    if (A.data()[i] != B.data()[i])
      return false;
  return true;
}

template <unsigned R, unsigned C> inline
bool operator!=( const MsqMatrix<R,C>& A, const MsqMatrix<R,C>& B ) 
  { return !(A == B); }

template <unsigned R, unsigned C> inline
std::ostream& operator<<( std::ostream& str, const MsqMatrix<R,C>& m )
{
  str << m.data()[0];
  for (unsigned i = 1; i < R*C; ++i)
    str << ' ' << m.data()[i];
  return str;
}

template <unsigned R, unsigned C> inline
std::istream& operator>>( std::istream& str, MsqMatrix<R,C>& m )
{
  for (unsigned i = 0; i < R*C; ++i)
    str >> m.data()[i];
  return str;
}

template <unsigned R> inline
double sqr_length( const MsqMatrix<R,1>& v )
  { return sqr_Frobenius(v); }

template <unsigned C> inline
double sqr_length( const MsqMatrix<1,C>& v )
  { return sqr_Frobenius(v); }

template <unsigned R> inline
double length( const MsqMatrix<R,1>& v )
  { return Frobenius(v); }

template <unsigned C> inline
double length( const MsqMatrix<1,C>& v )
  { return Frobenius(v); }

template <unsigned R, unsigned C> inline
double inner_product( const MsqMatrix<R,C>& m1, const MsqMatrix<R,C>& m2 )
{
  double result = 0.0;
  for (unsigned r = 0; r < R; ++r)
    for (unsigned c = 0; c < C; ++c)
      result += m1(r,c) * m2(r,c);
  return result;
}

template <unsigned R, unsigned C> inline
MsqMatrix<R,C> outer( const MsqMatrix<R,1>& v1, const MsqMatrix<C,1>& v2 )
{
  MsqMatrix<R,C> result;
  for (unsigned r = 0; r < R; ++r)
    for (unsigned c = 0; c < C; ++c)
      result(r,c) = v1(r,0) * v2(c,0);
  return result;
}

template <unsigned R, unsigned C> inline
MsqMatrix<R,C> outer( const MsqMatrix<1,R>& v1, const MsqMatrix<1,C>& v2 )
{
  MsqMatrix<R,C> result;
  for (unsigned r = 0; r < R; ++r)
    for (unsigned c = 0; c < C; ++c)
      result(r,c) = v1(0,r) * v2(0,c);
  return result;
}

inline
double vector_product( const MsqMatrix<2,1>& v1, const MsqMatrix<2,1>& v2 )
  { return v1(0,0) * v2(1,0) - v1(1,0) * v2(0,0); }

inline
double vector_product( const MsqMatrix<1,2>& v1, const MsqMatrix<1,2>& v2 )
  { return v1(0,0) * v2(0,1) - v1(0,1) * v2(0,0); }

inline
MsqMatrix<3,1> vector_product( const MsqMatrix<3,1>& a, const MsqMatrix<3,1>& b )
{
  MsqMatrix<3,1> result;
  result(0,0) = a(1,0)*b(2,0) - a(2,0)*b(1,0);
  result(1,0) = a(2,0)*b(0,0) - a(0,0)*b(2,0);
  result(2,0) = a(0,0)*b(1,0) - a(1,0)*b(0,0);
  return result;
}

inline
MsqMatrix<1,3> vector_product( const MsqMatrix<1,3>& a, const MsqMatrix<1,3>& b )
{
  MsqMatrix<1,3> result;
  result(0,0) = a(0,1)*b(0,2) - a(0,2)*b(0,1);
  result(0,1) = a(0,2)*b(0,0) - a(0,0)*b(0,2);
  result(0,2) = a(0,0)*b(0,1) - a(0,1)*b(0,0);
  return result;
}

template <unsigned R, unsigned C> inline
double operator%( const MsqMatrix<R,C>& v1, const MsqMatrix<R,C>& v2 )
  { return inner_product( v1, v2 ); }

inline
double operator*( const MsqMatrix<2,1>& v1, const MsqMatrix<2,1>& v2 )
  { return vector_product( v1, v2 ); }

inline
double operator*( const MsqMatrix<1,2>& v1, const MsqMatrix<1,2>& v2 )
  { return vector_product( v1, v2 ); }

inline
MsqMatrix<3,1> operator*( const MsqMatrix<3,1>& v1, const MsqMatrix<3,1>& v2 )
  { return vector_product( v1, v2 ); }

inline
MsqMatrix<1,3> operator*( const MsqMatrix<1,3>& v1, const MsqMatrix<1,3>& v2 )
  { return vector_product( v1, v2 ); }

/** Compute QR factorization of A */
inline
void QR( const MsqMatrix<3,3>& A, MsqMatrix<3,3>& Q, MsqMatrix<3,3>& R )
{
    Q = A;

    R(0,0) = sqrt(Q(0,0)*Q(0,0) + Q(1,0)*Q(1,0) + Q(2,0)*Q(2,0));
    double temp_dbl = 1.0/R(0,0);
    R(1,0) = 0.0;
    R(2,0) = 0.0;
    Q(0,0) *= temp_dbl;
    Q(1,0) *= temp_dbl;
    Q(2,0) *= temp_dbl;
    

    R(0,1)  = Q(0,0)*Q(0,1) + Q(1,0)*Q(1,1) + Q(2,0)*Q(2,1);
    Q(0,1) -= Q(0,0)*R(0,1);
    Q(1,1) -= Q(1,0)*R(0,1);
    Q(2,1) -= Q(2,0)*R(0,1);

    R(0,2)  = Q(0,0)*Q(0,2) + Q(1,0)*Q(1,2) + Q(2,0)*Q(2,2);
    Q(0,2) -= Q(0,0)*R(0,2);
    Q(1,2) -= Q(1,0)*R(0,2);
    Q(2,2) -= Q(2,0)*R(0,2);

    R(1,1) = sqrt(Q(0,1)*Q(0,1) + Q(1,1)*Q(1,1) + Q(2,1)*Q(2,1));
    temp_dbl = 1.0 / R(1,1);
    R(2,1) = 0.0;
    Q(0,1) *= temp_dbl;
    Q(1,1) *= temp_dbl;
    Q(2,1) *= temp_dbl;

    
    R(1,2)  = Q(0,1)*Q(0,2) + Q(1,1)*Q(1,2) + Q(2,1)*Q(2,2);
    Q(0,2) -= Q(0,1)*R(1,2);
    Q(1,2) -= Q(1,1)*R(1,2);
    Q(2,2) -= Q(2,1)*R(1,2);
  
    R(2,2) = sqrt(Q(0,2)*Q(0,2) + Q(1,2)*Q(1,2) + Q(2,2)*Q(2,2));
    temp_dbl = 1.0 / R(2,2);
    Q(0,2) *= temp_dbl;
    Q(1,2) *= temp_dbl;
    Q(2,2) *= temp_dbl;
}


/** Compute QR factorization of A */
inline
void QR( const MsqMatrix<2,2>& A, MsqMatrix<2,2>& Q, MsqMatrix<2,2>& R )
{
    R(0,0) = std::sqrt( A(0,0)*A(0,0) + A(1,0)*A(1,0) );
    const double r0inv = 1.0 / R(0,0);
    Q(0,0) = Q(1,1) = A(0,0) * r0inv;
    Q(1,0) = A(1,0) * r0inv;
    Q(0,1) = -Q(1,0);
    R(0,1) = r0inv * (A(0,0)*A(0,1) + A(1,0)*A(1,1));
    R(1,1) = r0inv * (A(0,0)*A(1,1) - A(0,1)*A(1,0));
    R(1,0) = 0.0;
}

} // namespace Mesquite

#endif
