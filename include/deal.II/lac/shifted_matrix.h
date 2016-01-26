// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2015 by the deal.II authors
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

#ifndef dealii__shifted_matrix_h
#define dealii__shifted_matrix_h


#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>

DEAL_II_NAMESPACE_OPEN

/*! @addtogroup Matrix2
 *@{
 */

/**
 * Matrix with shifted diagonal values.
 *
 * Given a matrix <tt>A</tt>, this class implements a matrix-vector product
 * with <i>A+s I</i>, where <i>s</i> is a provided shift parameter.
 *
 * @deprecated If deal.II was configured with C++11 support, use the
 * LinearOperator class instead, see the module on
 * @ref LAOperators "linear operators"
 * for further details.
 *
 * @author Guido Kanschat, 2000, 2001
 */
template<typename MatrixType>
class ShiftedMatrix
{
public:
  /**
   * Constructor.  Provide the base matrix and a shift parameter.
   */
  ShiftedMatrix (const MatrixType &A, const double sigma);

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
  template <typename VectorType>
  void vmult (VectorType &dst, const VectorType &src) const;

  /**
   * Residual.
   */
  template <typename VectorType>
  double residual (VectorType &dst, const VectorType &src, const VectorType &rhs) const;

private:
  /**
   * Storage for base matrix.
   */
  SmartPointer<const MatrixType,ShiftedMatrix<MatrixType> > A;

  /**
   * Auxiliary vector.
   */
  //    VectorType aux;
  /**
   * Shift parameter.
   */
  double sigma;
};



/**
 * Matrix with shifted diagonal values with respect to a certain scalar
 * product.
 *
 * Given a matrix <tt>A</tt>, this class implements a matrix-vector product
 * with <i>A+s M</i>, where <i>s</i> is a provided shift parameter and
 * <tt>M</tt> is the matrix representing the identity
 *
 * @deprecated If deal.II was configured with C++11 support, use the
 * LinearOperator class instead, see the module on
 * @ref LAOperators "linear operators"
 * for further details.
 *
 * @author Guido Kanschat, 2001
 */
template<typename MatrixType, class MatrixType2, class VectorType>
class ShiftedMatrixGeneralized
{
public:
  /**
   * Constructor. Provide the base matrix and a shift parameter.
   */
  ShiftedMatrixGeneralized (const MatrixType  &A,
                            const MatrixType2 &M,
                            const double       sigma);

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
  void vmult (VectorType &dst, const VectorType &src) const;

  /**
   * Residual.
   */
  double residual (VectorType &dst, const VectorType &src, const VectorType &rhs) const;

private:
  /**
   * Storage for base matrix.
   */
  SmartPointer<const MatrixType,ShiftedMatrixGeneralized<MatrixType,MatrixType2,VectorType> > A;
  /**
   * Storage for mass matrix.
   */
  SmartPointer<const MatrixType2,ShiftedMatrixGeneralized<MatrixType,MatrixType2,VectorType> > M;

  /**
   * Auxiliary vector.
   */
  mutable VectorType aux;

  /**
   * Shift parameter.
   */
  double sigma;
};


/*@}*/
//---------------------------------------------------------------------------

template <typename MatrixType>
inline
ShiftedMatrix<MatrixType>::ShiftedMatrix (const MatrixType &A,
                                          const double      sigma)
  :
  A(&A), sigma(sigma)
{}



template <typename MatrixType>
inline void
ShiftedMatrix<MatrixType>::shift (const double s)
{
  sigma = s;
}


template <typename MatrixType>
inline double
ShiftedMatrix<MatrixType>::shift () const
{
  return sigma;
}



template <typename MatrixType>
template <class VectorType>
inline void
ShiftedMatrix<MatrixType>::vmult (VectorType &dst, const VectorType &src) const
{
  A->vmult(dst, src);
  if (sigma != 0.)
    dst.add(sigma, src);
}


template <typename MatrixType>
template <class VectorType>
inline double
ShiftedMatrix<MatrixType>::residual (VectorType       &dst,
                                     const VectorType &src,
                                     const VectorType &rhs) const
{
  A->vmult(dst, src);
  if (sigma != 0.)
    dst.add(sigma, src);
  dst.sadd(-1.,1.,rhs);
  return dst.l2_norm ();
}


//---------------------------------------------------------------------------
template <typename MatrixType, class MatrixType2, class VectorType>
inline
ShiftedMatrixGeneralized<MatrixType, MatrixType2, VectorType>
::ShiftedMatrixGeneralized (const MatrixType  &A,
                            const MatrixType2 &M,
                            const double       sigma)
  :
  A(&A), M(&M), sigma(sigma)
{}


template <typename MatrixType, class MatrixType2, class VectorType>
inline void
ShiftedMatrixGeneralized<MatrixType, MatrixType2, VectorType>::shift (const double s)
{
  sigma = s;
}

template <typename MatrixType, class MatrixType2, class VectorType>
inline double
ShiftedMatrixGeneralized<MatrixType, MatrixType2, VectorType>::shift () const
{
  return sigma;
}


template <typename MatrixType, class MatrixType2, class VectorType>
inline void
ShiftedMatrixGeneralized<MatrixType, MatrixType2, VectorType>::vmult
(VectorType       &dst,
 const VectorType &src) const
{
  A->vmult(dst, src);
  if (sigma != 0.)
    {
      aux.reinit(dst);
      M->vmult(aux, src);
      dst.add(sigma, aux);
    }
}


template <typename MatrixType, class MatrixType2, class VectorType>
inline double
ShiftedMatrixGeneralized<MatrixType, MatrixType2, VectorType>::residual
(VectorType       &dst,
 const VectorType &src,
 const VectorType &rhs) const
{
  A->vmult(dst, src);
  if (sigma != 0.)
    {
      aux.reinit(dst);
      M->vmult(aux, src);
      dst.add(sigma, aux);
    }
  dst.sadd(-1.,1.,rhs);
  return dst.l2_norm ();
}


DEAL_II_NAMESPACE_CLOSE

#endif
