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

    size_type
    n_nonzero_elements() const;

    value_type
    el(const size_type i, const size_type j) const;

    /**
     * Set the element (i,j) to 'value'.
     */
    void
    set(const size_type i, const size_type j, const value_type value);

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

    void
    vmult(Vector &dst, const Vector &src) const;

    void
    Tvmult(Vector &dst, const Vector &src) const;

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

    MPI_Comm
    get_mpi_communicator() const;

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