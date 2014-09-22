// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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

#ifndef __deal2__transpose_matrix_h
#define __deal2__transpose_matrix_h


#include <deal.II/base/subscriptor.h>
#include <deal.II/lac/pointer_matrix.h>

DEAL_II_NAMESPACE_OPEN


/**
 * The transpose of a given matrix.  This auxiliary class swaps the
 * effect ov vmult() and Tvmult() as well as vmult_add() and
 * Tvmult_add().
 *
 * The implementation is analogous to the class PointerMatrix.
 *
 * @note The transposed matrix is never actually assembled. Instead,
 * only the matrix vector multiplication is performed in a transposed
 * way.
 *
 * @ingroup Matrix2
 * @author Guido Kanschat, 2006
 */
template<class MATRIX, class VECTOR>
class
  TransposeMatrix : public PointerMatrixBase<VECTOR>
{
public:
  /**
   * Constructor.  The pointer in the
   * argument is stored in this
   * class. As usual, the lifetime of
   * <tt>*M</tt> must be longer than the
   * one of the PointerMatrix.
   *
   * If <tt>M</tt> is zero, no
   * matrix is stored.
   */
  TransposeMatrix (const MATRIX *M=0);

  /**
   * Constructor. The name argument
   * is used to identify the
   * SmartPointer for this object.
   */
  TransposeMatrix(const char *name);

  /**
   * Constructor. <tt>M</tt> points
   * to a matrix which must live
   * longer than the
   * TransposeMatrix. The name
   * argument is used to identify
   * the SmartPointer for this
   * object.
   */
  TransposeMatrix(const MATRIX *M,
                  const char *name);

  // Use doc from base class
  virtual void clear();

  /**
   * Return whether the object is
   * empty.
   */
  bool empty () const;

  /**
   * Assign a new matrix
   * pointer. Deletes the old pointer
   * and releases its matrix.
   * @see SmartPointer
   */
  const TransposeMatrix &operator= (const MATRIX *M);

  /**
   * Matrix-vector product.
   */
  virtual void vmult (VECTOR &dst,
                      const VECTOR &src) const;

  /**
   * Tranposed matrix-vector product.
   */
  virtual void Tvmult (VECTOR &dst,
                       const VECTOR &src) const;

  /**
   * Matrix-vector product, adding to
   * <tt>dst</tt>.
   */
  virtual void vmult_add (VECTOR &dst,
                          const VECTOR &src) const;

  /**
   * Tranposed matrix-vector product,
   * adding to <tt>dst</tt>.
   */
  virtual void Tvmult_add (VECTOR &dst,
                           const VECTOR &src) const;

private:
  /**
   * The pointer to the actual matrix.
   */
  SmartPointer<const MATRIX,TransposeMatrix<MATRIX,VECTOR> > m;
};


//----------------------------------------------------------------------//


template<class MATRIX, class VECTOR>
TransposeMatrix<MATRIX, VECTOR>::TransposeMatrix (const MATRIX *M)
  : m(M)
{}


template<class MATRIX, class VECTOR>
TransposeMatrix<MATRIX, VECTOR>::TransposeMatrix (const char *name)
  : m(0, name)
{}


template<class MATRIX, class VECTOR>
TransposeMatrix<MATRIX, VECTOR>::TransposeMatrix (
  const MATRIX *M,
  const char *name)
  : m(M, name)
{}


template<class MATRIX, class VECTOR>
inline void
TransposeMatrix<MATRIX, VECTOR>::clear ()
{
  m = 0;
}


template<class MATRIX, class VECTOR>
inline const TransposeMatrix<MATRIX, VECTOR> &
TransposeMatrix<MATRIX, VECTOR>::operator= (const MATRIX *M)
{
  m = M;
  return *this;
}


template<class MATRIX, class VECTOR>
inline bool
TransposeMatrix<MATRIX, VECTOR>::empty () const
{
  if (m == 0)
    return true;
  return m->empty();
}

template<class MATRIX, class VECTOR>
inline void
TransposeMatrix<MATRIX, VECTOR>::vmult (VECTOR &dst,
                                        const VECTOR &src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->Tvmult (dst, src);
}


template<class MATRIX, class VECTOR>
inline void
TransposeMatrix<MATRIX, VECTOR>::Tvmult (VECTOR &dst,
                                         const VECTOR &src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->vmult (dst, src);
}


template<class MATRIX, class VECTOR>
inline void
TransposeMatrix<MATRIX, VECTOR>::vmult_add (VECTOR &dst,
                                            const VECTOR &src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->Tvmult_add (dst, src);
}


template<class MATRIX, class VECTOR>
inline void
TransposeMatrix<MATRIX, VECTOR>::Tvmult_add (VECTOR &dst,
                                             const VECTOR &src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->vmult_add (dst, src);
}


DEAL_II_NAMESPACE_CLOSE

#endif
