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

#ifndef dealii_psblas_sparsity_h
#define dealii_psblas_sparsity_h

#include <deal.II/lac/sparsity_pattern_base.h>

#ifdef DEAL_II_WITH_PSBLAS

#  include <deal.II/lac/psblas_common.h>

DEAL_II_NAMESPACE_OPEN

namespace PSCToolkitWrappers
{

  /**
   * Implementation of a distributed sparsity pattern based on PSBLAS
   * (Parallel Sparse BLAS), which is the
   * computational kernel of the
   * <a href="https://psctoolkit.github.io/">PSCToolkit</a> library.
   *
   * This class describes the non-zero structure of a distributed sparse
   * matrix.  Rows are partitioned across MPI processes according to an
   * IndexSet.  After construction the pattern is in the <em>Build</em>
   * state; calling compress() transitions it to the <em>Assembled</em>
   * state so that it can be passed to a SparseMatrix for allocation.
   *
   * @ingroup PSCToolkitWrappers
   * @ingroup Sparsity
   */
  class SparsityPattern : public SparsityPatternBase
  {
  public:
    using size_type = dealii::types::global_dof_index;

    /**
     * Default constructor.
     */
    SparsityPattern() = default;

    /**
     * Constructor of a square sparsity pattern from an existing IndexSet.
     */
    SparsityPattern(const IndexSet &parallel_partitioning,
                    const MPI_Comm  communicator = MPI_COMM_WORLD);

    /**
     * Destructor.
     */
    virtual ~SparsityPattern() override = default;

    /**
     * Add several elements in one row to the sparsity pattern.
     */
    template <typename ForwardIterator>
    void
    add_entries(const size_type row,
                ForwardIterator begin,
                ForwardIterator end,
                const bool      indices_are_sorted = false);

    /**
     * Add new entries to a given row. It calls the add_entries() function using
     * columns.begin() and columns.end() as begin and end iterators,
     * respectively.
     */
    virtual void
    add_row_entries(const size_type                  &row,
                    const ArrayView<const size_type> &columns,
                    const bool indices_are_sorted = false) override;

    /**
     * Add a single entry to the sparsity pattern.
     */
    void
    add(const size_type i, const size_type j);

    using SparsityPatternBase::add_entries;

    /**
     * This function compresses the sparsity pattern and allows the resulting
     * pattern to be used for actually generating a PSBLAS Sparse matrix.
     * This function must therefore be called once the structure is fixed.
     * Internally, it finalizes the descriptor object. This is a collective
     * operation, i.e., it needs to be run on all processors when used in
     * parallel.
     */
    void
    compress();

  private:
    /**
     * Shared pointer to the PSBLAS descriptor object.
     */
    std::shared_ptr<psb_c_descriptor> psblas_descriptor;

    friend class SparseMatrix;
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
