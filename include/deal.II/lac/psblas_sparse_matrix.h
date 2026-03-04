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

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparsity_pattern_base.h>
#include <deal.II/lac/vector.h>

#include <utility>


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
     * Construtor using an IndexSet and a MPI communicator to describe the
     * parallel partitioning of the matrix.
     */
    void
    reinit(const IndexSet &parallel_partitioning,
           const MPI_Comm  communicator = MPI_COMM_WORLD);

    /**
     * Construtor using an IndexSet and a MPI communicator to describe the
     * parallel partitioning of the matrix.
     */
    void
    reinit(const SparsityPattern &psblas_sparsity_pattern,
           const MPI_Comm         communicator = MPI_COMM_WORLD);

    size_type
    m() const;

    size_type
    n() const;

    size_type
    local_size() const;

    std::pair<size_type, size_type>
    local_range() const;

    bool
    in_local_range(const size_type index) const;

    size_type
    n_nonzero_elements() const;

    value_type
    el(const size_type i, const size_type j) const;

    value_type
    diag_element(const size_type i) const;

    /**
     * Set the element (i,j) to 'value'.
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
     * Add value 'value' to the element (i,j).
     */
    void
    add(const size_type i, const size_type j, const value_type value);

    void
    add(const size_type                row,
        const std::vector<size_type>  &col_indices,
        const std::vector<value_type> &values,
        const bool = false);

    void
    add(const size_type               row,
        const size_type               ncols,
        const std::vector<size_type> &col_indices,
        const value_type             *values,
        const bool = false,
        const bool = false);

    void
    add(const size_type                row,
        const size_type                ncols,
        const std::vector<size_type>  &col_indices,
        const std::vector<value_type> &values,
        const bool = false,
        const bool = false);

    void
    add(const size_type   row,
        const size_type   n_cols,
        const size_type  *col_indices,
        const value_type *values,
        const bool        elide_zero_values      = true,
        const bool        col_indices_are_sorted = false);

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
     * Get the underlying PSBLAS sparse matrix.
     */
    psb_c_dspmat *
    get_psblas_matrix() const;

    /**
     * Get the underlying PSBLAS descriptor.
     */
    psb_c_descriptor *
    get_psblas_descriptor() const;

    /**
     * Return the underlying MPI communicator.
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
    MPI_Comm communicator;

    psb_c_dspmat *psblas_sparse_matrix;

    std::shared_ptr<psb_c_descriptor> psblas_descriptor;

    psb_c_ctxt *psblas_context;
  };


} // namespace PSCToolkitWrappers


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PSBLAS
#endif // dealii_psctoolkit_h