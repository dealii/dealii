// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_psblas_matrix_h
#define dealii_psblas_matrix_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>


#ifdef DEAL_II_WITH_PSBLAS

#  include <deal.II/lac/psblas_sparsity_pattern.h>
#  include <deal.II/lac/psblas_vector.h>

DEAL_II_NAMESPACE_OPEN

namespace PSCToolkitWrappers
{

  class SparseMatrix : public EnableObserverPointer
  {
  public:
    /**
     * Type for container size.
     */
    using size_type = dealii::types::global_dof_index;

    /**
     * Type for container values.
     */
    using value_type = double;

    /**
     *Default constructor. Generates an empty (zero-size) matrix.
     */
    SparseMatrix();

    /**
     * Destructor. Internally, its frees the PSBLAS sparse matrix.
     */
    ~SparseMatrix();

    /**
     * Generate a matrix from a PSBLAS SparsityPattern.
     */
    SparseMatrix(const SparsityPattern &psblas_sparsity_pattern,
                 const MPI_Comm         communicator = MPI_COMM_WORLD);

    /**
     * Copy-constructor is deleted.
     *
     */
    SparseMatrix(const SparseMatrix &) = delete;

    /**
     * Copy assignment is deleted.
     *
     */
    SparseMatrix &
    operator=(const SparseMatrix &) = delete;

    /**
     * Assignment of a scalar to all elements of the matrix is deleted.
     */
    SparseMatrix &
    operator=(const value_type d) = delete;

    /**
     * Copy the contents of another SparseMatrix into this one. Currently not
     * implemented from PSBLAS.
     */
    void
    copy_from(const SparseMatrix &other);


    /**
     * Initialize using an IndexSet and a MPI communicator to describe the
     * parallel partitioning of the matrix.
     */
    void
    reinit(const IndexSet &parallel_partitioning,
           const MPI_Comm  communicator = MPI_COMM_WORLD);

    /**
     * Initialize using an IndexSet and a MPI communicator to describe the
     * parallel partitioning of the matrix.
     */
    void
    reinit(const SparsityPattern &psblas_sparsity_pattern,
           const MPI_Comm         communicator = MPI_COMM_WORLD);

    /**
     * Initialize a square matrix where the size() of the IndexSet determines
     * the number of rows and columns.
     */
    void
    reinit(const IndexSet               &local_rows,
           const DynamicSparsityPattern &sparsity_pattern,
           const MPI_Comm                communicator = MPI_COMM_WORLD);

    /**
     * Return the number of rows in this matrix.
     */
    size_type
    m() const;

    /**
     * Return the number of columns in this matrix.
     */
    size_type
    n() const;

    /**
     * Return the local dimension of the matrix, i.e. the number of rows stored
     * on the present MPI process. For sequential matrices, this number is the
     * same as m(), but for parallel matrices it may be smaller.
     * To figure out which elements exactly are stored locally, use
     * local_range().
     */
    size_type
    local_size() const;

    /**
     * Return a pair of indices indicating which rows of this matrix are
     * stored locally. The first number is the index of the first row stored,
     * the second the index of the one past the last one that is stored
     * locally. If this is a sequential matrix, then the result will be the
     * pair (0,m()), otherwise it will be a pair (i,i+n), where
     * <tt>n=local_size()</tt>.
     */
    std::pair<size_type, size_type>
    local_range() const;

    /**
     * Return whether @p index is in the local range or not, see also
     * local_range().
     */
    bool
    in_local_range(const size_type index) const;

    /**
     * Return the number of nonzero elements of this matrix.
     */
    size_type
    n_nonzero_elements() const;

    /**
     * Return the value of the matrix entry (<i>i,j</i>).
     *
     * @note Currently not implemented for PSBLAS.
     */
    value_type
    el(const size_type i, const size_type j) const;

    /**
     * Return the diagonal element of the matrix at row <i>i</i>.
     *
     * @note Currently not implemented for PSBLAS.
     */
    value_type
    diag_element(const size_type i) const;

    /**
     * Return the element (i,j) of the matrix. This function is equivalent to
     * the el() function.
     *
     * @note Currently not implemented for PSBLAS.
     */
    value_type
    operator()(const size_type i, const size_type j) const;

    /**
     * Set the element (i,j) to 'value'. If <tt>value</tt> is
     * not a finite number an exception is thrown.
     */
    void
    set(const size_type i, const size_type j, const value_type value);

    /**
     * Set all elements given in a FullMatrix into the sparse matrix locations
     * given by <tt>indices</tt>. In other words, this function writes the
     * elements in <tt>full_matrix</tt> into the calling matrix, using the
     * local-to-global indexing specified by <tt>indices</tt> for both the rows
     * and the columns of the matrix. This function assumes a quadratic sparse
     * matrix and a quadratic full_matrix, the usual situation in FE
     * calculations.
     */
    void
    set(const std::vector<size_type> &indices,
        const FullMatrix<double>     &matrix);

    /**
     * Add value 'value' to the element (i,j). If <tt>value</tt> is
     * not a finite number an exception is thrown.
     */
    void
    add(const size_type i, const size_type j, const value_type value);

    /**
     * Set several elements in the specified row of the matrix with column
     * indices as given by @p col_indices to the respective value.
     */
    void
    add(const size_type                row,
        const std::vector<size_type>  &col_indices,
        const std::vector<value_type> &values,
        const bool = false);

    /**
     * Add an array of values given by @p values in the given global matrix @p row at
     * columns specified by @p col_indices in the sparse matrix.
     */
    void
    add(const size_type               row,
        const size_type               ncols,
        const std::vector<size_type> &col_indices,
        const value_type             *values,
        const bool = false,
        const bool = false);

    /**
     * Same as above, but with a vector of values instead of a pointer to an
     * array. The vector must have the same size as the number of columns
     given by @p ncols and the size of @p col_indices.
     */
    void
    add(const size_type                row,
        const size_type                ncols,
        const std::vector<size_type>  &col_indices,
        const std::vector<value_type> &values,
        const bool = false,
        const bool = false);

    /**
     * Same as above, but using raw pointers for the column indices and values.
     * The array of column indices must have the same size as @p ncols, and
     * the array of values must have the same size as @p ncols as well.
     */
    void
    add(const size_type   row,
        const size_type   n_cols,
        const size_type  *col_indices,
        const value_type *values,
        const bool        elide_zero_values      = true,
        const bool        col_indices_are_sorted = false);

    /**
     * Compress the matrix after all insertions have been done.
     */
    void
    compress();

    /**
     * Matrix-vector multiplication: let <i>dst = M*src</i> with <i>M</i>
     * being this matrix.
     *
     * Source and destination must not be the same vector.
     */
    void
    vmult(Vector &dst, const Vector &src) const;

    /**
     * Adding matrix-vector multiplication: Add <i>M*src</i> to <i>dst</i> with
     * <i>M</i> being this matrix.
     *
     * Source and destination must not be the same vector.
     */
    void
    vmult_add(Vector &dst, const Vector &src) const;


    /**
     * Matrix-vector multiplication: let <i>dst = M<sup>T</sup>*src</i> with
     * <i>M</i> being this matrix. This function does the same as vmult() but
     * takes the transposed matrix.
     *
     * Source and destination must not be the same vector.
     */
    void
    Tvmult(Vector &dst, const Vector &src) const;

    /**
     * Adding matrix-vector multiplication: Add <i>M^T*src</i> to <i>dst</i>
     * with <i>M</i> being this matrix.
     *
     * Source and destination must not be the same vector.
     */
    void
    Tvmult_add(Vector &dst, const Vector &src) const;

    /**
     * Get the underlying PSBLAS sparse matrix. Do not use this function unless
     * you know what you are doing.
     */
    psb_c_dspmat *
    get_psblas_matrix() const;

    /**
     * Get the underlying PSBLAS descriptor. Do not use this function unless you
     * know what you are doing.
     */
    psb_c_descriptor *
    get_psblas_descriptor() const;

    /**
     * Return the underlying MPI communicator. Do not use this function unless
     * you know what you are doing.
     */
    MPI_Comm
    get_mpi_communicator() const;


    /**
     * Exception
     */
    DeclExceptionMsg(ExcSourceEqualsDestination,
                     "You are attempting an operation on two vectors that "
                     "are the same object, but the operation requires that the "
                     "two objects are in fact different.");

  private:
    /**
     * The underlying MPI communicator.
     *
     */
    MPI_Comm communicator;

    /**
     * Pointer the underlying PSBLAS sparse matrix.
     *
     */
    psb_c_dspmat *psblas_sparse_matrix;

    /**
     * Shared pointer to the underlying PSBLAS descriptor object.
     *
     */
    std::shared_ptr<psb_c_descriptor> psblas_descriptor;

    /**
     * Pointer to the underlying PSBLAS context object.
     *
     */
    psb_c_ctxt *psblas_context;

    /**
     * State of the descriptor associated with the vector. Its state can be
     * either default, building or assembled).
     */
    internal::State state;

    friend class PreconditionAMG;
  };


} // namespace PSCToolkitWrappers


DEAL_II_NAMESPACE_CLOSE
#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PSBLAS
#endif
