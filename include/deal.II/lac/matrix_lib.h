// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_matrix_lib_h
#define dealii_matrix_lib_h

#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/vector_memory.h>

DEAL_II_NAMESPACE_OPEN

template <typename number>
class Vector;
template <typename number>
class BlockVector;
template <typename number>
class SparseMatrix;

/*! @addtogroup Matrix2
 *@{
 */


/**
 * Mean value filter.  The vmult() functions of this matrix filter out mean
 * values of the vector.  If the vector is of type BlockVector, then an
 * additional parameter selects a single component for this operation.
 *
 * In mathematical terms, this class acts as if it was the matrix $I-\frac
 * 1n{\mathbf 1}_n{\mathbf 1}_n^T$ where ${\mathbf 1}_n$ is a vector of size
 * $n$ that has only ones as its entries. Thus, taking the dot product between
 * a vector $\mathbf v$ and $\frac 1n {\mathbf 1}_n$ yields the <i>mean
 * value</i> of the entries of ${\mathbf v}$. Consequently, $ \left[I-\frac
 * 1n{\mathbf 1}_n{\mathbf 1}_n^T\right] \mathbf v = \mathbf v - \left[\frac
 * 1n {\mathbf v} \cdot {\mathbf 1}_n\right]{\mathbf 1}_n$ subtracts from every
 * vector element the mean value of all elements.
 *
 * @deprecated Use a LinearOperator, or a BlockLinearOperator instead. you
 * can construct such a filter by using mean_value_filter, or
 * block_diagonal_operator (with a mean_value_filter block), respectively.
 *
 * @author Guido Kanschat, 2002, 2003
 */
class DEAL_II_DEPRECATED MeanValueFilter : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Constructor, optionally selecting a component.
   */
  MeanValueFilter(const size_type component = numbers::invalid_size_type);

  /**
   * Subtract mean value from @p v.
   */
  template <typename number>
  void
  filter(Vector<number> &v) const;

  /**
   * Subtract mean value from @p v.
   */
  template <typename number>
  void
  filter(BlockVector<number> &v) const;

  /**
   * Return the source vector with subtracted mean value.
   */
  template <typename number>
  void
  vmult(Vector<number> &dst, const Vector<number> &src) const;

  /**
   * Add source vector with subtracted mean value to dest.
   */
  template <typename number>
  void
  vmult_add(Vector<number> &dst, const Vector<number> &src) const;

  /**
   * Return the source vector with subtracted mean value in selected
   * component.
   */
  template <typename number>
  void
  vmult(BlockVector<number> &dst, const BlockVector<number> &src) const;

  /**
   * Add a source to dest, where the mean value in the selected component is
   * subtracted.
   */
  template <typename number>
  void
  vmult_add(BlockVector<number> &dst, const BlockVector<number> &src) const;


  /**
   * Not implemented.
   */
  template <typename VectorType>
  void
  Tvmult(VectorType &, const VectorType &) const;

  /**
   * Not implemented.
   */
  template <typename VectorType>
  void
  Tvmult_add(VectorType &, const VectorType &) const;

private:
  /**
   * Component for filtering block vectors.
   */
  const size_type component;
};



/*@}*/
//---------------------------------------------------------------------------


template <typename VectorType>
inline void
MeanValueFilter::Tvmult(VectorType &, const VectorType &) const
{
  Assert(false, ExcNotImplemented());
}


template <typename VectorType>
inline void
MeanValueFilter::Tvmult_add(VectorType &, const VectorType &) const
{
  Assert(false, ExcNotImplemented());
}


DEAL_II_NAMESPACE_CLOSE

#endif
