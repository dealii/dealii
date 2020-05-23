// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2019 by the deal.II authors
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

#ifndef dealii_diagonal_matrix_h
#define dealii_diagonal_matrix_h


#include <deal.II/base/config.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_operation.h>

DEAL_II_NAMESPACE_OPEN

/**
 * This class represents a <i>n x n</i> diagonal matrix based on a vector of
 * size <i>n</i>. The matrix-vector products are realized by @p
 * VectorType::scale, so the template vector class needs to provide a
 * @p scale() method.
 *
 * When using this class with ConstraintsMatrix::distribute_local_to_global(),
 * the underlying vector needs to provide write access to all entries referenced
 * by cells in an assembly process. This means that this class also needs access
 * to ghost entries that are owned by other processors than the calling one.
 * In practice this requires initialization of the vector as follows
 * @code
 * DiagonalMatrix<LinearAlgebra::distributed::Vector<double> > diagonal_matrix;
 * LinearAlgebra::distributed::Vector<double> &diagonal_vector =
 *   diagonal_matrix.get_vector();
 * diagonal_vector.reinit(locally_owned_dofs,
 *                        locally_relevant_dofs,
 *                        mpi_communicator);
 * @endcode
 *
 * @author Martin Kronbichler, 2016
 */
template <typename VectorType = Vector<double>>
class DiagonalMatrix : public Subscriptor
{
public:
  using value_type = typename VectorType::value_type;
  using size_type  = typename VectorType::size_type;

  /**
   * Default constructor. The object needs still to be reinitialized to be
   * usable.
   */
  DiagonalMatrix() = default;

  /**
   * Constructor initializing this object as a diagonal matrix of size `n x n`
   * where `n` is the size of the vector, and with diagonal entries equal to the
   * elements of @p vec.
   */
  explicit DiagonalMatrix(const VectorType &vec);

  /**
   * Initialize with a given vector by copying the content of the vector
   * @p vec.
   */
  void
  reinit(const VectorType &vec);

  /**
   * Compresses the data structures and allows the resulting matrix to be used
   * in all other operations like matrix-vector products. This is a collective
   * operation, i.e., it needs to be run on all processors when used in
   * parallel.
   */
  void
  compress(VectorOperation::values operation);

  /**
   * Return a reference to the underlying vector for manipulation of the
   * entries on the matrix diagonal.
   */
  VectorType &
  get_vector();

  /**
   * Clear content of this object and reset to the state of default constructor.
   */
  void
  clear();

  /**
   * Return a read-only reference to the underlying vector.
   */
  const VectorType &
  get_vector() const;

  /**
   * Number of rows of this matrix. This number corresponds to the size of the
   * underlying vector.
   */
  size_type
  m() const;

  /**
   * Number of columns of this matrix. This number corresponds to the size of
   * the underlying vector.
   */
  size_type
  n() const;

  /**
   * Read-only access to a value. This is restricted to the case where
   * <i>i==j</i> due to the matrix storage.
   *
   * If the vector representing the diagonal is distributed with MPI, not all
   * of the indices <i>i</i> might actually be accessible. Refer to the method
   * <code>get_vector().locally_owned_elements()</code> for the entries that
   * actually are accessible.
   */
  value_type
  operator()(const size_type i, const size_type j) const;

  /**
   * Read-write access to a value. This is restricted to the case where
   * <i>i==j</i> due to the matrix storage.
   *
   * If the vector representing the diagonal is distributed with MPI, not all
   * of the indices <i>i</i> might actually be accessible. Refer to the method
   * <code>get_vector().locally_owned_elements()</code> for the entries that
   * actually are accessible.
   */
  value_type &
  operator()(const size_type i, const size_type j);

  /**
   * Add an array of values given by <tt>values</tt> in the given global
   * matrix row at columns specified by col_indices. Due to the storage of
   * this matrix, entries are only added to the diagonal of the matrix. All
   * other entries are ignored and no exception is thrown.
   *
   * This function is for a consistent interface with the other matrix
   * classes in deal.II and can be used in
   * AffineConstraints::distribute_local_to_global to get exactly the same
   * diagonal as when assembling into a sparse matrix.
   */
  template <typename number2>
  void
  add(const size_type  row,
      const size_type  n_cols,
      const size_type *col_indices,
      const number2 *  values,
      const bool       elide_zero_values      = true,
      const bool       col_indices_are_sorted = false);

  /**
   * Add value to the element (i,j).
   *
   * Due to the storage of this matrix, entries are only added to the diagonal
   * of the matrix. All other entries are ignored and no exception is thrown.
   */
  void
  add(const size_type i, const size_type j, const value_type value);

  /**
   * Performs a matrix-vector multiplication with the given matrix.
   */
  void
  vmult(VectorType &dst, const VectorType &src) const;

  /**
   * Performs a transpose matrix-vector multiplication with the given
   * matrix. Since this represents a diagonal matrix, exactly the same as
   * vmult().
   */
  void
  Tvmult(VectorType &dst, const VectorType &src) const;

  /**
   * Adds the result of a matrix-vector multiplication into the destination
   * vector dst. Needs to create a temporary vector, which makes performance
   * slower than for @p vmult().
   */
  void
  vmult_add(VectorType &dst, const VectorType &src) const;

  /**
   * Adds the result of a transpose matrix-vector multiplication into the
   * destination vector dst. Needs to create a temporary vector, which makes
   * performance slower than for @p Tvmult().
   */
  void
  Tvmult_add(VectorType &dst, const VectorType &src) const;

  /**
   * Initialize vector @p dst to have the same size and partition as
   * @p diagonal member of this class.
   *
   * This is a part of the interface required
   * by linear_operator().
   */
  void
  initialize_dof_vector(VectorType &dst) const;

  /**
   * Return the memory consumption of this object.
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * The stored vector.
   */
  VectorType diagonal;
};

/* ---------------------------------- Inline functions ------------------- */

#ifndef DOXYGEN

template <typename VectorType>
DiagonalMatrix<VectorType>::DiagonalMatrix(const VectorType &vec)
  : diagonal(vec)
{}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::clear()
{
  diagonal.reinit(0);
}



template <typename VectorType>
std::size_t
DiagonalMatrix<VectorType>::memory_consumption() const
{
  return diagonal.memory_consumption();
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::reinit(const VectorType &vec)
{
  diagonal = vec;
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::initialize_dof_vector(VectorType &dst) const
{
  dst.reinit(diagonal);
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::compress(VectorOperation::values operation)
{
  diagonal.compress(operation);
}



template <typename VectorType>
VectorType &
DiagonalMatrix<VectorType>::get_vector()
{
  return diagonal;
}



template <typename VectorType>
const VectorType &
DiagonalMatrix<VectorType>::get_vector() const
{
  return diagonal;
}



template <typename VectorType>
typename VectorType::size_type
DiagonalMatrix<VectorType>::m() const
{
  return diagonal.size();
}



template <typename VectorType>
typename VectorType::size_type
DiagonalMatrix<VectorType>::n() const
{
  return diagonal.size();
}



template <typename VectorType>
typename VectorType::value_type
DiagonalMatrix<VectorType>::operator()(const size_type i,
                                       const size_type j) const
{
  Assert(i == j, ExcIndexRange(j, i, i + 1));
  (void)j;
  return diagonal(i);
}



template <typename VectorType>
typename VectorType::value_type &
DiagonalMatrix<VectorType>::operator()(const size_type i, const size_type j)
{
  Assert(i == j, ExcIndexRange(j, i, i + 1));
  (void)j;
  return diagonal(i);
}



template <typename VectorType>
template <typename number2>
void
DiagonalMatrix<VectorType>::add(const size_type  row,
                                const size_type  n_cols,
                                const size_type *col_indices,
                                const number2 *  values,
                                const bool,
                                const bool)
{
  for (size_type i = 0; i < n_cols; ++i)
    if (col_indices[i] == row)
      diagonal(row) += values[i];
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::add(const size_type  i,
                                const size_type  j,
                                const value_type value)
{
  if (i == j)
    diagonal(i) += value;
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::vmult(VectorType &dst, const VectorType &src) const
{
  dst = src;
  dst.scale(diagonal);
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::Tvmult(VectorType &dst, const VectorType &src) const
{
  vmult(dst, src);
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::vmult_add(VectorType &      dst,
                                      const VectorType &src) const
{
  VectorType tmp(src);
  tmp.scale(diagonal);
  dst += tmp;
}



template <typename VectorType>
void
DiagonalMatrix<VectorType>::Tvmult_add(VectorType &      dst,
                                       const VectorType &src) const
{
  vmult_add(dst, src);
}


#endif

DEAL_II_NAMESPACE_CLOSE

#endif
