// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#ifndef __deal2__matrix_iterator_h
#define __deal2__matrix_iterator_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

DEAL_II_NAMESPACE_OPEN

/**
 * STL conforming iterator for constant
 * and non-constant matrices.
 *
 * This iterator is abstracted from
 * the actual matrix type and can be
 * used for any matrix having the
 * required ACCESSOR type.
 *
 * @author Guido Kanschat, 2006, based on previous a implementation
 */
template <class ACCESSOR>
class MatrixIterator
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Typedef for the matrix type
   * (including constness) we are to
   * operate on.
   */
  typedef typename ACCESSOR::MatrixType MatrixType;

  /**
   * Constructor. Create an
   * iterator into the matrix
   * <tt>matrix</tt> for the given
   * <tt>row</tt> and the
   * <tt>index</tt> within it.
   */
  MatrixIterator (MatrixType      *matrix,
                  const size_type  row = 0,
                  const size_type  index = 0);

  /**
   * Copy from another matrix
   * iterator. Mostly implemented
   * to allow initialization of a
   * constant iterator from a non
   * constant, this function only
   * requires that a conversion
   * from the other iterator's
   * accessor to this accessor
   * object is possible.
   */
  template <class OtherAccessor>
  MatrixIterator(const MatrixIterator<OtherAccessor> &other);

  /**
   * Prefix increment.
   */
  MatrixIterator &operator++ ();

  /**
   * Postfix increment.
   */
  MatrixIterator operator++ (int);

  /**
   * Dereferencing operator.
   */
  const ACCESSOR &operator* () const;

  /**
   * Dereferencing operator.
   */
  const ACCESSOR *operator-> () const;

  /**
   * Comparison. True, if
   * both accessors are equal.
   */
  bool operator == (const MatrixIterator &) const;

  /**
   * Inverse of <tt>==</tt>.
   */
  bool operator != (const MatrixIterator &) const;

  /**
   * Comparison operator. Result is
   * true if either the first row
   * number is smaller or if the row
   * numbers are equal and the first
   * index is smaller.
   *
   * This function is only valid if
   * both iterators point into the same
   * matrix.
   */
  bool operator < (const MatrixIterator &) const;

  /**
   * Comparison operator. Works in the
   * same way as above operator, just
   * the other way round.
   */
  bool operator > (const MatrixIterator &) const;

private:
  /**
   * Store an object of the
   * accessor class.
   */
  ACCESSOR accessor;

  /**
   * Allow other iterators access
   * to private data.
   */
  template <class OtherAccessor> friend class MatrixIterator;
};


//----------------------------------------------------------------------//

template <class ACCESSOR>
inline
MatrixIterator<ACCESSOR>::
MatrixIterator (MatrixType      *matrix,
                const size_type  r,
                const size_type  i)
  :
  accessor(matrix, r, i)
{}


template <class ACCESSOR>
template <class OtherAccessor>
inline
MatrixIterator<ACCESSOR>::
MatrixIterator (const MatrixIterator<OtherAccessor> &other)
  :
  accessor(other.accessor)
{}


template <class ACCESSOR>
inline
MatrixIterator<ACCESSOR> &
MatrixIterator<ACCESSOR>::operator++ ()
{
  accessor.advance ();
  return *this;
}


template <class ACCESSOR>
inline
MatrixIterator<ACCESSOR>
MatrixIterator<ACCESSOR>::operator++ (int)
{
  const MatrixIterator iter = *this;
  accessor.advance ();
  return iter;
}


template <class ACCESSOR>
inline
const ACCESSOR &
MatrixIterator<ACCESSOR>::operator* () const
{
  return accessor;
}


template <class ACCESSOR>
inline
const ACCESSOR *
MatrixIterator<ACCESSOR>::operator-> () const
{
  return &accessor;
}


template <class ACCESSOR>
inline
bool
MatrixIterator<ACCESSOR>::
operator == (const MatrixIterator &other) const
{
  return (accessor == other.accessor);
}


template <class ACCESSOR>
inline
bool
MatrixIterator<ACCESSOR>::
operator != (const MatrixIterator &other) const
{
  return ! (*this == other);
}


template <class ACCESSOR>
inline
bool
MatrixIterator<ACCESSOR>::
operator < (const MatrixIterator &other) const
{
  Assert (&accessor.get_matrix() == &other.accessor.get_matrix(),
          ExcInternalError());

  return (accessor < other.accessor);
}


template <class ACCESSOR>
inline
bool
MatrixIterator<ACCESSOR>::
operator > (const MatrixIterator &other) const
{
  return (other < *this);
}

DEAL_II_NAMESPACE_CLOSE

#endif
