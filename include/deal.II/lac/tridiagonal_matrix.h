// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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

#ifndef __deal2__tridiagonal_matrix_h
#define __deal2__tridiagonal_matrix_h

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/lac/lapack_support.h>

#include <vector>
#include <iomanip>

DEAL_II_NAMESPACE_OPEN

// forward declarations
template<typename number> class Vector;


/*! @addtogroup Matrix1
 *@{
 */


/**
 * A quadratic tridiagonal matrix. That is, a matrix where all entries
 * are zero, except the diagonal and the entries left and right of it.
 *
 * The matrix has an additional symmetric mode, in which case only the
 * upper triangle of the matrix is stored and mirrored to the lower
 * one for matrix vector operations.
 *
 * @ingroup Matrix1
 * @author Guido Kanschat, 2005, 2006
*/
template<typename number>
class TridiagonalMatrix
{
public:
  ///@name Constructors
  //@{
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * @name Constructors and initalization.
   */
  /**
   * Constructor generating an
   * empty matrix of dimension
   * <tt>n</tt>.
   */
  TridiagonalMatrix(size_type    n = 0,
                    bool symmetric = false);

  /**
   * Reinitialize the matrix to a
   * new size and reset all entries
   * to zero. The symmetry
   * properties may be set as well.
   */
  void reinit(size_type n,
              bool symmetric = false);


//@}
///@name Non-modifying operators
//@{

  /**
   * Number of rows of this matrix.
   * To remember: this matrix is an
   * <i>m x m</i>-matrix.
   */
  size_type m () const;

  /**
   * Number of columns of this matrix.
   * To remember: this matrix is an
   * <i>n x n</i>-matrix.
   */
  size_type n () const;

  /**
   * Return whether the matrix
   * contains only elements with
   * value zero. This function is
   * mainly for internal
   * consistency checks and should
   * seldom be used when not in
   * debug mode since it uses quite
   * some time.
   */
  bool all_zero () const;



//@}
///@name Element access
//@{
  /**
   * Read-only access to a
   * value. This is restricted to
   * the case where <i>|i-j| <=
   * 1</i>.
   */
  number operator()(size_type i, size_type j) const;

  /**
   * Read-write access to a
   * value. This is restricted to
   * the case where <i>|i-j| <=
   * 1</i>.
   *
   * @note In case of symmetric
   * storage technique, the entries
   * <i>(i,j)</i> and <i>(j,i)</i>
   * are identified and <b>both</b>
   * exist. This must be taken into
   * account if adding up is used
   * for matrix assembling in order
   * not to obtain doubled entries.
   */
  number &operator()(size_type i, size_type j);

//@}
///@name Multiplications with vectors
//@{

  /**
   * Matrix-vector-multiplication. Multiplies
   * <tt>v</tt> from the right and
   * stores the result in
   * <tt>w</tt>.
   *
   * If the optional parameter
   * <tt>adding</tt> is <tt>true</tt>, the
   * result is added to <tt>w</tt>.
   *
   * Source and destination must
   * not be the same vector.
   */
  void vmult (Vector<number>       &w,
              const Vector<number> &v,
              const bool            adding=false) const;

  /**
   * Adding
   * Matrix-vector-multiplication. Same
   * as vmult() with parameter
   * <tt>adding=true</tt>, but
   * widely used in
   * <tt>deal.II</tt> classes.
   *
   * Source and destination must
   * not be the same vector.
   */
  void vmult_add (Vector<number>       &w,
                  const Vector<number> &v) const;

  /**
   * Transpose
   * matrix-vector-multiplication.
   * Multiplies
   * <tt>v<sup>T</sup></tt> from
   * the left and stores the result
   * in <tt>w</tt>.
   *
   * If the optional parameter
   * <tt>adding</tt> is <tt>true</tt>, the
   * result is added to <tt>w</tt>.
   *
   * Source and destination must
   * not be the same vector.
   */
  void Tvmult (Vector<number>       &w,
               const Vector<number> &v,
               const bool            adding=false) const;

  /**
   * Adding transpose
   * matrix-vector-multiplication. Same
   * as Tvmult() with parameter
   * <tt>adding=true</tt>, but
   * widely used in
   * <tt>deal.II</tt> classes.
   *
   * Source and destination must
   * not be the same vector.
   */
  void Tvmult_add (Vector<number>       &w,
                   const Vector<number> &v) const;

  /**
   * Build the matrix scalar product
   * <tt>u^T M v</tt>. This function is mostly
   * useful when building the cellwise
   * scalar product of two functions in
   * the finite element context.
   */
  number matrix_scalar_product (const Vector<number> &u,
                                const Vector<number> &v) const;

  /**
   * Return the square of the norm
   * of the vector <tt>v</tt> with
   * respect to the norm induced by
   * this matrix,
   * i.e. <i>(v,Mv)</i>. This is
   * useful, e.g. in the finite
   * element context, where the
   * <i>L<sup>2</sup></i> norm of a
   * function equals the matrix
   * norm with respect to the mass
   * matrix of the vector
   * representing the nodal values
   * of the finite element
   * function.
   *
   * Obviously, the matrix needs to
   * be quadratic for this operation.
   */
  number matrix_norm_square (const Vector<number> &v) const;

//@}
///@name Matrixnorms
//@{

  /**
   * Return the $l_1$-norm of the matrix, i.e.
   * $|M|_1=max_{all columns j}\sum_{all
   * rows i} |M_ij|$,
   * (max. sum of columns). This is the
   * natural matrix norm that is compatible
   * to the $l_1$-norm for vectors, i.e.
   * $|Mv|_1\leq |M|_1 |v|_1$.
   * (cf. Rannacher Numerik0)
   */
  number l1_norm () const;

  /**
   * Return the $l_\infty$-norm of the
   * matrix, i.e.
   * $|M|_\infty=\max_{all rows i}\sum_{all
   * columns j} |M_{ij}|$,
   * (max. sum of rows).
   * This is the
   * natural matrix norm that is compatible
   * to the $l_\infty$-norm of vectors, i.e.
   * $|Mv|_\infty \leq |M|_\infty |v|_\infty$.
   */
  number linfty_norm () const;

  /**
   * The Frobenius norm of the matrix.
   * Return value is the root of the square
   * sum of all matrix entries.
   */
  number frobenius_norm () const;

  /**
   * Compute the relative norm of
   * the skew-symmetric part. The
   * return value is the Frobenius
   * norm of the skew-symmetric
   * part of the matrix divided by
   * that of the matrix.
   *
   * Main purpose of this function
   * is to check, if a matrix is
   * symmetric within a certain
   * accuracy, or not.
   */
  number relative_symmetry_norm2 () const;
//@}
///@name LAPACK operations
//@{
  /**
   * Compute the eigenvalues of the
   * symmetric tridiagonal matrix.
   *
   * @note This function requires
   * configuration of deal.II with
   * LAPACK support. Additionally,
   * the matrix must use symmetric
   * storage technique.
   */
  void compute_eigenvalues();
  /**
   * After calling
   * compute_eigenvalues(), you can
   * access each eigenvalue here.
   */
  number eigenvalue(const size_type i) const;
//@}
///@name Miscellanea
//@{
  /**
   * Output of the matrix in
   * user-defined format.
   */
  template <class OUT>
  void print (OUT &s,
              const unsigned int  width=5,
              const unsigned int  precision=2) const;

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t memory_consumption () const;
  //@}

private:
  /**
   * The diagonal entries.
   */
  std::vector<number> diagonal;
  /**
   * The entries left of the
   * diagonal. The entry with index
   * zero is always zero, since the
   * first row has no entry left of
   * the diagonal. Therefore, the
   * length of this vector is the
   * same as that of #diagonal.
   *
   * The length of this vector is
   * zero for symmetric storage. In
   * this case, the second element
   * of #left is identified with
   * the first element of #right.
   */
  std::vector<number> left;
  /**
   * The entries right of the
   * diagonal. The last entry is
   * always zero, since the last
   * row has no entry right of the
   * diagonal. Therefore, the
   * length of this vector is the
   * same as that of #diagonal.
   */
  std::vector<number> right;

  /**
   * If this flag is true, only the
   * entries to the right of the
   * diagonal are stored and the
   * matrix is assumed symmetric.
   */
  bool is_symmetric;

  /**
   * The state of the
   * matrix. Normally, the state is
   * LAPACKSupport::matrix,
   * indicating that the object can
   * be used for regular matrix
   * operations.
   *
   * See explanation of this data
   * type for details.
   */
  LAPACKSupport::State state;
};

/**@}*/

//---------------------------------------------------------------------------
#ifndef DOXYGEN

template<typename number>
types::global_dof_index
TridiagonalMatrix<number>::m() const
{
  return diagonal.size();
}



template<typename number>
types::global_dof_index
TridiagonalMatrix<number>::n() const
{
  return diagonal.size();
}


template<typename number>
inline
number
TridiagonalMatrix<number>::operator()(size_type i, size_type j) const
{
  Assert(i<n(), ExcIndexRange(i,0,n()));
  Assert(j<n(), ExcIndexRange(j,0,n()));
  Assert (i<=j+1, ExcIndexRange(i,j-1,j+2));
  Assert (j<=i+1, ExcIndexRange(j,i-1,i+2));

  if (j==i)
    return diagonal[i];
  if (j==i-1)
    {
      if (is_symmetric)
        return right[i-1];
      else
        return left[i];
    }

  if (j==i+1)
    return right[i];

  Assert (false, ExcInternalError());
  return 0;
}


template<typename number>
inline
number &
TridiagonalMatrix<number>::operator()(size_type i, size_type j)
{
  Assert(i<n(), ExcIndexRange(i,0,n()));
  Assert(j<n(), ExcIndexRange(j,0,n()));
  Assert (i<=j+1, ExcIndexRange(i,j-1,j+2));
  Assert (j<=i+1, ExcIndexRange(j,i-1,i+2));

  if (j==i)
    return diagonal[i];
  if (j==i-1)
    {
      if (is_symmetric)
        return right[i-1];
      else
        return left[i];
    }

  if (j==i+1)
    return right[i];

  Assert (false, ExcInternalError());
  return diagonal[0];
}


template <typename number>
template <class OUT>
void
TridiagonalMatrix<number>::print (
  OUT &s,
  const unsigned int width,
  const unsigned int) const
{
  for (size_type i=0; i<n(); ++i)
    {
      if (i>0)
        s << std::setw(width) << (*this)(i,i-1);
      else
        s << std::setw(width) << "";

      s << ' ' << (*this)(i,i) << ' ';

      if (i<n()-1)
        s << std::setw(width) << (*this)(i,i+1);

      s << std::endl;
    }
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

