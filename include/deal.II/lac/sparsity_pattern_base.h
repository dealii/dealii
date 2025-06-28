// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_sparsity_pattern_base_h
#define dealii_sparsity_pattern_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/types.h>

#include <utility>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Sparsity
 * @{
 */

/**
 * Base class for all sparsity patterns, defining a common interface by which
 * new values can be added.
 */
class SparsityPatternBase : public EnableObserverPointer
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Constructor. Sets up an empty (zero-by-zero) sparsity pattern.
   */
  SparsityPatternBase();

  /**
   * Constructor. Sets up a @p rows by @p cols sparsity pattern.
   */
  SparsityPatternBase(const size_type rows, const size_type cols);

  /**
   * Copy constructor.
   */
  SparsityPatternBase(const SparsityPatternBase &sparsity_pattern) = default;

  /**
   * Move constructor.
   */
  SparsityPatternBase(SparsityPatternBase &&sparsity_pattern) noexcept =
    default;

  /**
   * Assignment operator.
   */
  SparsityPatternBase &
  operator=(const SparsityPatternBase &sparsity_pattern) = default;

  /**
   * Move assignment operator.
   */
  SparsityPatternBase &
  operator=(SparsityPatternBase &&sparsity_pattern) noexcept = default;

  /**
   * Return number of rows of this matrix, which equals the dimension of the
   * image space.
   */
  size_type
  n_rows() const;

  /**
   * Return number of columns of this matrix, which equals the dimension of
   * the range space.
   */
  size_type
  n_cols() const;

  /**
   * Optimized function for adding new entries to a given row.
   */
  virtual void
  add_row_entries(const size_type                  &row,
                  const ArrayView<const size_type> &columns,
                  const bool indices_are_sorted = false) = 0;

  /**
   * General function for adding new entries from an unsorted list of pairs.
   *
   * This function is useful when multiple entries need to be added which do
   * not correspond to a particular row, e.g., when assembling a flux sparsity
   * pattern.
   */
  virtual void
  add_entries(const ArrayView<const std::pair<size_type, size_type>> &entries);

protected:
  /**
   * Internal function for updating the stored size of the sparsity pattern.
   */
  virtual void
  resize(const size_type rows, const size_type cols);

  /**
   * Number of rows that this sparsity pattern shall represent.
   */
  size_type rows;

  /**
   * Number of columns that this sparsity pattern shall represent.
   */
  size_type cols;
};

/**
 * @}
 */


/* ---------------------------- Inline functions ---------------------------- */

#ifndef DOXYGEN

inline SparsityPatternBase::SparsityPatternBase()
  : rows(0)
  , cols(0)
{}



inline SparsityPatternBase::SparsityPatternBase(const size_type rows,
                                                const size_type cols)
  : rows(rows)
  , cols(cols)
{}



inline SparsityPatternBase::size_type
SparsityPatternBase::n_rows() const
{
  return rows;
}



inline SparsityPatternBase::size_type
SparsityPatternBase::n_cols() const
{
  return cols;
}



inline void
SparsityPatternBase::resize(const size_type n_rows, const size_type n_cols)
{
  this->rows = n_rows;
  this->cols = n_cols;
}
#endif

DEAL_II_NAMESPACE_CLOSE

#endif
