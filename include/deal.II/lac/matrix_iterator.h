// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_matrix_iterator_h
#define dealii_matrix_iterator_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/types.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Iterator for constant and non-constant matrices.
 *
 * This iterator is abstracted from the actual matrix type and can be used for
 * any matrix having the required AccessorType type.
 */
template <typename AccessorType>
class MatrixIterator
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Typedef for the matrix type (including constness) we are to operate on.
   */
  using MatrixType = typename AccessorType::MatrixType;

  /**
   * Constructor. Create an iterator into the matrix <tt>matrix</tt> for the
   * given <tt>row</tt> and the <tt>index</tt> within it.
   */
  MatrixIterator(MatrixType     *matrix,
                 const size_type row   = 0,
                 const size_type index = 0);

  /**
   * Copy from another matrix iterator. Mostly implemented to allow
   * initialization of a constant iterator from a non constant, this function
   * only requires that a conversion from the other iterator's AccessorType to
   * this AccessorType object is possible.
   */
  template <class OtherAccessorType>
  MatrixIterator(const MatrixIterator<OtherAccessorType> &other);

  /**
   * Prefix increment.
   */
  MatrixIterator &
  operator++();

  /**
   * Postfix increment.
   */
  MatrixIterator
  operator++(int);

  /**
   * Dereferencing operator.
   */
  const AccessorType &
  operator*() const;

  /**
   * Dereferencing operator.
   */
  const AccessorType *
  operator->() const;

  /**
   * Comparison. True, if both accessors are equal.
   */
  template <class OtherAccessorType>
  bool
  operator==(const MatrixIterator<OtherAccessorType> &) const;

  /**
   * Inverse of <tt>==</tt>.
   */
  template <class OtherAccessorType>
  bool
  operator!=(const MatrixIterator<OtherAccessorType> &) const;

  /**
   * Comparison operator. Result is true if either the first row number is
   * smaller or if the row numbers are equal and the first index is smaller.
   *
   * This function is only valid if both iterators point into the same matrix.
   */
  bool
  operator<(const MatrixIterator &) const;

  /**
   * Comparison operator. Works in the same way as above operator, just the
   * other way round.
   */
  bool
  operator>(const MatrixIterator &) const;

private:
  /**
   * Store an object of the AccessorType class.
   */
  AccessorType accessor;

  // Allow other iterators access to private data.
  template <class OtherAccessorType>
  friend class MatrixIterator;
};


//----------------------------------------------------------------------//

template <typename AccessorType>
inline MatrixIterator<AccessorType>::MatrixIterator(MatrixType     *matrix,
                                                    const size_type r,
                                                    const size_type i)
  : accessor(matrix, r, i)
{}


template <typename AccessorType>
template <class OtherAccessorType>
inline MatrixIterator<AccessorType>::MatrixIterator(
  const MatrixIterator<OtherAccessorType> &other)
  : accessor(other.accessor)
{}


template <typename AccessorType>
inline MatrixIterator<AccessorType> &
MatrixIterator<AccessorType>::operator++()
{
  accessor.advance();
  return *this;
}


template <typename AccessorType>
inline MatrixIterator<AccessorType>
MatrixIterator<AccessorType>::operator++(int)
{
  const MatrixIterator iter = *this;
  accessor.advance();
  return iter;
}


template <typename AccessorType>
inline const AccessorType &
MatrixIterator<AccessorType>::operator*() const
{
  return accessor;
}


template <typename AccessorType>
inline const AccessorType *
MatrixIterator<AccessorType>::operator->() const
{
  return &accessor;
}


template <typename AccessorType>
template <class OtherAccessorType>
inline bool
MatrixIterator<AccessorType>::operator==(
  const MatrixIterator<OtherAccessorType> &other) const
{
  return (accessor == other.accessor);
}


template <typename AccessorType>
template <class OtherAccessorType>
inline bool
MatrixIterator<AccessorType>::operator!=(
  const MatrixIterator<OtherAccessorType> &other) const
{
  return !(*this == other);
}


template <typename AccessorType>
inline bool
MatrixIterator<AccessorType>::operator<(const MatrixIterator &other) const
{
  Assert(&accessor.get_matrix() == &other.accessor.get_matrix(),
         ExcInternalError());

  return (accessor < other.accessor);
}


template <typename AccessorType>
inline bool
MatrixIterator<AccessorType>::operator>(const MatrixIterator &other) const
{
  return (other < *this);
}

DEAL_II_NAMESPACE_CLOSE

#endif
