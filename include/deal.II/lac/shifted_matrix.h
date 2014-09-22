// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__shifted_matrix_h
#define __deal2__shifted_matrix_h


#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>

DEAL_II_NAMESPACE_OPEN

/*! @addtogroup Matrix2
 *@{
 */

/**
 * Matrix with shifted diagonal values.
 *
 * Given a matrix <tt>A</tt>, this class implements a matrix-vector
 * product with <i>A+s I</i>, where <i>s</i> is a provided shift
 * parameter.
 *
 * @author Guido Kanschat, 2000, 2001
 */
template<class MATRIX>
class ShiftedMatrix
{
public:
  /**
   * Constructor.  Provide the base
   * matrix and a shift parameter.
   */
  ShiftedMatrix (const MATRIX &A, const double sigma);

  /**
   * Set the shift parameter.
   */
  void shift (const double sigma);

  /**
   * Access to the shift parameter.
   */
  double shift () const;

  /**
   * Matrix-vector-product.
   */
  template <class VECTOR>
  void vmult (VECTOR &dst, const VECTOR &src) const;

  /**
   * Residual.
   */
  template <class VECTOR>
  double residual (VECTOR &dst, const VECTOR &src, const VECTOR &rhs) const;

private:
  /**
   * Storage for base matrix.
   */
  SmartPointer<const MATRIX,ShiftedMatrix<MATRIX> > A;

  /**
   * Auxiliary vector.
   */
  //    VECTOR aux;
  /**
   * Shift parameter.
   */
  double sigma;
};



/**
 * Matrix with shifted diagonal values with respect to a certain scalar product.
 *
 * Given a matrix <tt>A</tt>, this class implements a matrix-vector product
 * with <i>A+s M</i>, where <i>s</i> is a provided shift parameter and
 * <tt>M</tt> is the matrix representing the identity
 *
 * @author Guido Kanschat, 2001
 */
template<class MATRIX, class MASSMATRIX, class VECTOR>
class ShiftedMatrixGeneralized
{
public:
  /**
   * Constructor.
   * Provide the base matrix and a shift parameter.
   */
  ShiftedMatrixGeneralized (const MATRIX &A,
                            const MASSMATRIX &M,
                            const double sigma);

  /**
   * Set the shift parameter.
   */
  void shift (const double sigma);

  /**
   * Access to the shift parameter.
   */
  double shift () const;

  /**
   * Matrix-vector-product.
   */
  void vmult (VECTOR &dst, const VECTOR &src) const;

  /**
   * Residual.
   */
  double residual (VECTOR &dst, const VECTOR &src, const VECTOR &rhs) const;

private:
  /**
   * Storage for base matrix.
   */
  SmartPointer<const MATRIX,ShiftedMatrixGeneralized<MATRIX,MASSMATRIX,VECTOR> > A;
  /**
   * Storage for mass matrix.
   */
  SmartPointer<const MASSMATRIX,ShiftedMatrixGeneralized<MATRIX,MASSMATRIX,VECTOR> > M;

  /**
   * Auxiliary vector.
   */
  VECTOR aux;

  /**
   * Shift parameter.
   */
  double sigma;
};


/*@}*/
//---------------------------------------------------------------------------

template <class MATRIX>
inline
ShiftedMatrix<MATRIX>::ShiftedMatrix (const MATRIX &A, const double sigma)
  :
  A(&A), sigma(sigma)
{}



template <class MATRIX>
inline void
ShiftedMatrix<MATRIX>::shift (const double s)
{
  sigma = s;
}


template <class MATRIX>
inline double
ShiftedMatrix<MATRIX>::shift () const
{
  return sigma;
}



template <class MATRIX>
template <class VECTOR>
inline void
ShiftedMatrix<MATRIX>::vmult (VECTOR &dst, const VECTOR &src) const
{
  A->vmult(dst, src);
  if (sigma != 0.)
    dst.add(sigma, src);
}


template <class MATRIX>
template <class VECTOR>
inline double
ShiftedMatrix<MATRIX>::residual (VECTOR &dst,
                                 const VECTOR &src,
                                 const VECTOR &rhs) const
{
  A->vmult(dst, src);
  if (sigma != 0.)
    dst.add(sigma, src);
  dst.sadd(-1.,1.,rhs);
  return dst.l2_norm ();
}


//---------------------------------------------------------------------------
template <class MATRIX, class MASSMATRIX, class VECTOR>
inline
ShiftedMatrixGeneralized<MATRIX, MASSMATRIX, VECTOR>
::ShiftedMatrixGeneralized (const MATRIX &A,
                            const MASSMATRIX &M,
                            const double sigma)
  :
  A(&A), M(&M), sigma(sigma)
{}


template <class MATRIX, class MASSMATRIX, class VECTOR>
inline void
ShiftedMatrixGeneralized<MATRIX, MASSMATRIX, VECTOR>::shift (const double s)
{
  sigma = s;
}

template <class MATRIX, class MASSMATRIX, class VECTOR>
inline double
ShiftedMatrixGeneralized<MATRIX, MASSMATRIX, VECTOR>::shift () const
{
  return sigma;
}


template <class MATRIX, class MASSMATRIX, class VECTOR>
inline void
ShiftedMatrixGeneralized<MATRIX, MASSMATRIX, VECTOR>::vmult (VECTOR &dst,
    const VECTOR &src) const
{
  A->vmult(dst, src);
  if (sigma != 0.)
    {
      aux.reinit(dst);
      M.vmult(aux, src);
      dst.add(sigma, aux);
    }
}


template <class MATRIX, class MASSMATRIX, class VECTOR>
inline double
ShiftedMatrixGeneralized<MATRIX, MASSMATRIX, VECTOR>::residual (VECTOR &dst,
    const VECTOR &src,
    const VECTOR &rhs) const
{
  A->vmult(dst, src);
  if (sigma != 0.)
    {
      aux.reinit(dst);
      M.vmult(aux, src);
      dst.add(sigma, aux);
    }
  dst.sadd(-1.,1.,rhs);
  return dst.l2_norm ();
}


DEAL_II_NAMESPACE_CLOSE

#endif
