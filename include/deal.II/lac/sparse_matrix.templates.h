// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_sparse_matrix_templates_h
#define dealii_sparse_matrix_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/parallel.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <boost/io/ios_state.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <numeric>
#include <ostream>
#include <vector>



DEAL_II_NAMESPACE_OPEN


template <typename number>
SparseMatrix<number>::SparseMatrix()
  : cols(nullptr, "SparseMatrix")
  , val(nullptr)
  , max_len(0)
{}



template <typename number>
SparseMatrix<number>::SparseMatrix(const SparseMatrix &m)
  : EnableObserverPointer(m)
  , cols(nullptr, "SparseMatrix")
  , val(nullptr)
  , max_len(0)
{
  Assert(m.cols == nullptr && m.val == nullptr && m.max_len == 0,
         ExcMessage(
           "This constructor can only be called if the provided argument "
           "is an empty matrix. This constructor can not be used to "
           "copy-construct a non-empty matrix. Use the "
           "SparseMatrix::copy_from() function for that purpose."));
}



template <typename number>
SparseMatrix<number>::SparseMatrix(SparseMatrix<number> &&m) noexcept
  : EnableObserverPointer(std::move(m))
  , cols(std::move(m.cols))
  , val(std::move(m.val))
  , max_len(m.max_len)
{
  m.cols    = nullptr;
  m.val     = nullptr;
  m.max_len = 0;
}



template <typename number>
SparseMatrix<number> &
SparseMatrix<number>::operator=(const SparseMatrix<number> &m)
{
  Assert(m.cols == nullptr && m.val == nullptr && m.max_len == 0,
         ExcMessage("This operator can only be called if the provided right "
                    "hand side is an empty matrix. This operator can not be "
                    "used to copy a non-empty matrix. Use the "
                    "SparseMatrix::copy_from() function for that purpose."));

  return *this;
}



template <typename number>
SparseMatrix<number> &
SparseMatrix<number>::operator=(SparseMatrix<number> &&m) noexcept
{
  cols    = m.cols;
  val     = std::move(m.val);
  max_len = m.max_len;

  m.cols    = nullptr;
  m.val     = nullptr;
  m.max_len = 0;

  return *this;
}



template <typename number>
SparseMatrix<number>::SparseMatrix(const SparsityPattern &c)
  : cols(nullptr, "SparseMatrix")
  , val(nullptr)
  , max_len(0)
{
  // virtual functions called in constructors and destructors never use the
  // override in a derived class
  // for clarity be explicit on which function is called
  SparseMatrix<number>::reinit(c);
}



template <typename number>
SparseMatrix<number>::SparseMatrix(const SparsityPattern &c,
                                   const IdentityMatrix  &id)
  : cols(nullptr, "SparseMatrix")
  , val(nullptr)
  , max_len(0)
{
  Assert(c.n_rows() == id.m(), ExcDimensionMismatch(c.n_rows(), id.m()));
  Assert(c.n_cols() == id.n(), ExcDimensionMismatch(c.n_cols(), id.n()));

  // virtual functions called in constructors and destructors never use the
  // override in a derived class
  // for clarity be explicit on which function is called
  SparseMatrix<number>::reinit(c);
  for (size_type i = 0; i < n(); ++i)
    this->set(i, i, 1.);
}



template <typename number>
SparseMatrix<number>::~SparseMatrix()
{
  cols = nullptr;
}



template <typename number>
SparseMatrix<number> &
SparseMatrix<number>::operator=(const double d)
{
  Assert(d == 0, ExcScalarAssignmentOnlyForZeroValue());

  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(cols->compressed || cols->empty(),
         SparsityPattern::ExcNotCompressed());

  // do initial zeroing of elements in parallel. Try to achieve a similar
  // layout as when doing matrix-vector products, as on some NUMA systems, a
  // memory block is assigned to memory banks where the first access is
  // generated. For sparse matrices, the first operations is usually the
  // operator=. The grain size is chosen to reflect the number of rows in
  // minimum_parallel_grain_size, weighted by the number of nonzero entries
  // per row on average.
  const std::size_t matrix_size = cols->n_nonzero_elements();
  const size_type   grain_size =
    internal::SparseMatrixImplementation::minimum_parallel_grain_size *
    (cols->n_nonzero_elements() + m()) / m();
  if (matrix_size > grain_size)
    parallel::apply_to_subranges(
      0U,
      matrix_size,
      [values = this->val.get()](const size_type begin, const size_type end) {
        std::fill(values + begin, values + end, number(0.));
      },
      grain_size);
  else if (matrix_size > 0)
    {
      std::fill(val.get(), val.get() + matrix_size, 0);
    }

  return *this;
}



template <typename number>
SparseMatrix<number> &
SparseMatrix<number>::operator=(const IdentityMatrix &id)
{
  Assert(cols->n_rows() == id.m(),
         ExcDimensionMismatch(cols->n_rows(), id.m()));
  Assert(cols->n_cols() == id.n(),
         ExcDimensionMismatch(cols->n_cols(), id.n()));

  *this = 0;
  for (size_type i = 0; i < n(); ++i)
    this->set(i, i, 1.);

  return *this;
}



template <typename number>
void
SparseMatrix<number>::reinit(const SparsityPattern &sparsity)
{
  cols = &sparsity;

  if (cols->empty())
    {
      val.reset();
      max_len = 0;
      return;
    }

  const std::size_t N = cols->n_nonzero_elements();
  if (N > max_len || max_len == 0)
    {
      val     = std::make_unique<number[]>(N);
      max_len = N;
    }

  *this = 0.;
}



template <typename number>
void
SparseMatrix<number>::clear()
{
  cols = nullptr;
  val.reset();
  max_len = 0;
}



template <typename number>
bool
SparseMatrix<number>::empty() const
{
  if (cols == nullptr)
    return true;
  else
    return cols->empty();
}



template <typename number>
typename SparseMatrix<number>::size_type
SparseMatrix<number>::get_row_length(const size_type row) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  return cols->row_length(row);
}



template <typename number>
std::size_t
SparseMatrix<number>::n_nonzero_elements() const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  return cols->n_nonzero_elements();
}



template <typename number>
std::size_t
SparseMatrix<number>::n_actually_nonzero_elements(const double threshold) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(threshold >= 0, ExcMessage("Negative threshold!"));
  size_type         nnz       = 0;
  const std::size_t nnz_alloc = n_nonzero_elements();
  for (size_type i = 0; i < nnz_alloc; ++i)
    if (std::abs(val[i]) > threshold)
      ++nnz;
  return nnz;
}



template <typename number>
void
SparseMatrix<number>::symmetrize()
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(cols->rows == cols->cols, ExcNotQuadratic());

  const size_type n_rows = m();
  for (size_type row = 0; row < n_rows; ++row)
    {
      // first skip diagonal entry
      number *val_ptr = &val[cols->rowstart[row]];
      if (m() == n())
        ++val_ptr;
      const size_type    *colnum_ptr = &cols->colnums[cols->rowstart[row] + 1];
      const number *const val_end_of_row = &val[cols->rowstart[row + 1]];

      // treat lower left triangle
      while ((val_ptr != val_end_of_row) && (*colnum_ptr < row))
        {
          // compute the mean of this
          // and the transpose value
          const number mean_value =
            (*val_ptr + val[(*cols)(*colnum_ptr, row)]) / number(2.0);
          // set this value and the
          // transpose one to the
          // mean
          *val_ptr = mean_value;
          set(*colnum_ptr, row, mean_value);

          // advance pointers
          ++val_ptr;
          ++colnum_ptr;
        };
    };
}



template <typename number>
template <typename somenumber>
SparseMatrix<number> &
SparseMatrix<number>::copy_from(const SparseMatrix<somenumber> &matrix)
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  Assert(cols == matrix.cols, ExcDifferentSparsityPatterns());

  std::copy(matrix.val.get(),
            matrix.val.get() + cols->n_nonzero_elements(),
            val.get());

  return *this;
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::copy_from(const FullMatrix<somenumber> &matrix)
{
  // first delete previous content
  *this = 0;

  // then copy old matrix
  for (size_type row = 0; row < matrix.m(); ++row)
    for (size_type col = 0; col < matrix.n(); ++col)
      if (matrix(row, col) != somenumber())
        set(row, col, number(matrix(row, col)));
}



#ifdef DEAL_II_WITH_TRILINOS

template <typename number>
SparseMatrix<number> &
SparseMatrix<number>::copy_from(const TrilinosWrappers::SparseMatrix &matrix)
{
  Assert(m() == matrix.m(), ExcDimensionMismatch(m(), matrix.m()));
  Assert(n() == matrix.n(), ExcDimensionMismatch(n(), matrix.n()));

  // first delete previous content
  *this = 0;

  std::vector<TrilinosScalar> value_cache;
  std::vector<size_type>      colnum_cache;

  for (size_type row = 0; row < matrix.m(); ++row)
    {
      value_cache.resize(matrix.n());
      colnum_cache.resize(matrix.n());

      // copy column indices and values and at the same time enquire about the
      // length of the row
      int ncols;
      int ierr = matrix.trilinos_matrix().ExtractGlobalRowCopy(
        row,
        matrix.row_length(row),
        ncols,
        value_cache.data(),
        reinterpret_cast<TrilinosWrappers::types::int_type *>(
          colnum_cache.data()));
      Assert(ierr == 0, ExcTrilinosError(ierr));

      // resize arrays to the size actually used
      value_cache.resize(ncols);
      colnum_cache.resize(ncols);

      // then copy everything in one swoop
      this->set(row, colnum_cache, value_cache);
    }

  return *this;
}

#endif


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::add(const number                    factor,
                          const SparseMatrix<somenumber> &matrix)
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  Assert(cols == matrix.cols, ExcDifferentSparsityPatterns());

  number             *val_ptr    = val.get();
  const somenumber   *matrix_ptr = matrix.val.get();
  const number *const end_ptr    = val.get() + cols->n_nonzero_elements();

  while (val_ptr != end_ptr)
    *val_ptr++ += factor * number(*matrix_ptr++);
}



namespace internal
{
  namespace SparseMatrixImplementation
  {
    using size_type = types::global_dof_index;

    /**
     * Perform a vmult using the SparseMatrix data structures, but only using
     * a subinterval for the row indices.
     *
     * In the sequential case, this function is called on all rows, in the
     * parallel case it may be called on a subrange, at the discretion of the
     * task scheduler.
     */
    template <typename number, typename InVector, typename OutVector>
    void
    vmult_on_subrange(const size_type    begin_row,
                      const size_type    end_row,
                      const number      *values,
                      const std::size_t *rowstart,
                      const size_type   *colnums,
                      const InVector    &src,
                      OutVector         &dst,
                      const bool         add)
    {
      const number                *val_ptr    = &values[rowstart[begin_row]];
      const size_type             *colnum_ptr = &colnums[rowstart[begin_row]];
      typename OutVector::iterator dst_ptr    = dst.begin() + begin_row;

      if (add == false)
        for (size_type row = begin_row; row < end_row; ++row)
          {
            typename OutVector::value_type s   = 0.;
            const number *const val_end_of_row = &values[rowstart[row + 1]];
            while (val_ptr != val_end_of_row)
              s += typename OutVector::value_type(*val_ptr++) *
                   typename OutVector::value_type(src(*colnum_ptr++));
            *dst_ptr++ = s;
          }
      else
        for (size_type row = begin_row; row < end_row; ++row)
          {
            typename OutVector::value_type s   = *dst_ptr;
            const number *const val_end_of_row = &values[rowstart[row + 1]];
            while (val_ptr != val_end_of_row)
              s += typename OutVector::value_type(*val_ptr++) *
                   typename OutVector::value_type(src(*colnum_ptr++));
            *dst_ptr++ = s;
          }
    }
  } // namespace SparseMatrixImplementation
} // namespace internal


template <typename number>
template <typename number2>
void
SparseMatrix<number>::add(const size_type  row,
                          const size_type  n_cols,
                          const size_type *col_indices,
                          const number2   *values,
                          const bool       elide_zero_values,
                          const bool       col_indices_are_sorted)
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  AssertIndexRange(row, m());

  // if we have sufficiently many columns
  // and sorted indices it is faster to
  // just go through the column indices and
  // look whether we found one, rather than
  // doing many binary searches
  if (elide_zero_values == false && col_indices_are_sorted == true &&
      n_cols > 3)
    {
      // check whether the given indices are
      // really sorted
      if constexpr (running_in_debug_mode())
        {
          for (size_type i = 1; i < n_cols; ++i)
            Assert(col_indices[i] > col_indices[i - 1],
                   ExcMessage(
                     "List of indices is unsorted or contains duplicates."));
        }

      const size_type *this_cols    = &cols->colnums[cols->rowstart[row]];
      const size_type  row_length_1 = cols->row_length(row) - 1;
      number          *val_ptr      = &val[cols->rowstart[row]];

      if (m() == n())
        {
          // find diagonal and add it if found
          Assert(this_cols[0] == row, ExcInternalError());
          const size_type *diag_pos =
            Utilities::lower_bound(col_indices, col_indices + n_cols, row);
          const size_type diag      = diag_pos - col_indices;
          size_type       post_diag = diag;
          if (diag != n_cols && *diag_pos == row)
            {
              val_ptr[0] += *(values + (diag_pos - col_indices));
              ++post_diag;
            }

          // Add indices before diagonal. Because the input array
          // is sorted, and because the entries in this matrix row
          // are sorted, we can just linearly walk the colnums array
          // and the input array in parallel, stopping whenever the
          // former matches the column index of the next index in
          // the input array:
          size_type counter = 1;
          for (size_type i = 0; i < diag; ++i)
            {
              while (this_cols[counter] < col_indices[i] &&
                     counter < row_length_1)
                ++counter;

              Assert((this_cols[counter] == col_indices[i]) ||
                       (values[i] == number2()),
                     ExcInvalidIndex(row, col_indices[i]));

              val_ptr[counter] += values[i];
            }

          // Then do the same to add indices after the diagonal:
          for (size_type i = post_diag; i < n_cols; ++i)
            {
              while (this_cols[counter] < col_indices[i] &&
                     counter < row_length_1)
                ++counter;

              Assert((this_cols[counter] == col_indices[i]) ||
                       (values[i] == number2()),
                     ExcInvalidIndex(row, col_indices[i]));

              val_ptr[counter] += values[i];
            }

          Assert(counter < cols->row_length(row),
                 ExcMessage(
                   "Specified invalid column indices in add function."));
        }
      else
        {
          // Use the same algorithm as above, but because the matrix is
          // not square, we can now do without the split for diagonal/
          // entries before the diagonal/entries are the diagonal.
          size_type counter = 0;
          for (size_type i = 0; i < n_cols; ++i)
            {
              while (this_cols[counter] < col_indices[i] &&
                     counter < row_length_1)
                ++counter;

              Assert((this_cols[counter] == col_indices[i]) ||
                       (values[i] == number2()),
                     ExcInvalidIndex(row, col_indices[i]));

              val_ptr[counter] += values[i];
            }
          Assert(counter < cols->row_length(row),
                 ExcMessage(
                   "Specified invalid column indices in add function."));
        }
      return;
    }

  // unsorted case: first, search all the
  // indices to find out which values we
  // actually need to add.
  const size_type *const my_cols        = cols->colnums.get();
  size_type              index          = cols->rowstart[row];
  const size_type        next_row_index = cols->rowstart[row + 1];

  for (size_type j = 0; j < n_cols; ++j)
    {
      const number value = number(values[j]);
      AssertIsFinite(value);

      if constexpr (running_in_debug_mode())
        {
          if (elide_zero_values == true && value == number())
            continue;
        }
      else
        {
          if (value == number())
            continue;
        }

      // check whether the next index to add is
      // the next present index in the sparsity
      // pattern (otherwise, do a binary
      // search)
      if (index < next_row_index && my_cols[index] == col_indices[j])
        goto add_value;

      index = cols->operator()(row, col_indices[j]);

      // it is allowed to add elements to
      // the matrix that are not part of
      // the sparsity pattern, if the
      // value we add is zero
      if (index == SparsityPattern::invalid_entry)
        {
          Assert(value == number(), ExcInvalidIndex(row, col_indices[j]));
          continue;
        }

    add_value:
      val[index] += value;
      ++index;
    }
}



template <typename number>
template <typename number2>
void
SparseMatrix<number>::set(const size_type  row,
                          const size_type  n_cols,
                          const size_type *col_indices,
                          const number2   *values,
                          const bool       elide_zero_values)
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  AssertIndexRange(row, m());

  // First, search all the indices to find
  // out which values we actually need to
  // set.
  const size_type  *my_cols = cols->colnums.get();
  std::size_t       index = cols->rowstart[row], next_index = index;
  const std::size_t next_row_index = cols->rowstart[row + 1];

  if (elide_zero_values == true)
    {
      for (size_type j = 0; j < n_cols; ++j)
        {
          const number value = number(values[j]);
          AssertIsFinite(value);

          if (value == number())
            continue;

          // check whether the next index to set is
          // the next present index in the sparsity
          // pattern (otherwise, do a binary
          // search)
          if (index != next_row_index && my_cols[index] == col_indices[j])
            goto set_value;

          next_index = cols->operator()(row, col_indices[j]);

          // it is allowed to set elements in
          // the matrix that are not part of
          // the sparsity pattern, if the
          // value to which we set it is zero
          if (next_index == SparsityPattern::invalid_entry)
            {
              Assert(false, ExcInvalidIndex(row, col_indices[j]));
              continue;
            }
          index = next_index;

        set_value:
          val[index] = value;
          ++index;
        }
    }
  else
    {
      // same code as above, but now check for zeros
      for (size_type j = 0; j < n_cols; ++j)
        {
          const number value = number(values[j]);
          AssertIsFinite(value);

          if (index != next_row_index && my_cols[index] == col_indices[j])
            goto set_value_checked;

          next_index = cols->operator()(row, col_indices[j]);

          if (next_index == SparsityPattern::invalid_entry)
            {
              Assert(value == number(), ExcInvalidIndex(row, col_indices[j]));
              continue;
            }
          index = next_index;

        set_value_checked:
          val[index] = value;
          ++index;
        }
    }
}



template <typename number>
template <class OutVector, class InVector>
void
SparseMatrix<number>::vmult(OutVector &dst, const InVector &src) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  Assert(m() == dst.size(), ExcDimensionMismatch(m(), dst.size()));
  Assert(n() == src.size(), ExcDimensionMismatch(n(), src.size()));

  Assert(!PointerComparison::equal(&src, &dst), ExcSourceEqualsDestination());

  parallel::apply_to_subranges(
    0U,
    m(),
    [this, &src, &dst](const size_type begin_row, const size_type end_row) {
      internal::SparseMatrixImplementation::vmult_on_subrange(
        begin_row,
        end_row,
        val.get(),
        cols->rowstart.get(),
        cols->colnums.get(),
        src,
        dst,
        false);
    },
    internal::SparseMatrixImplementation::minimum_parallel_grain_size);
}



template <typename number>
template <class OutVector, class InVector>
void
SparseMatrix<number>::Tvmult(OutVector &dst, const InVector &src) const
{
  Assert(val != nullptr, ExcNotInitialized());
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(n() == dst.size(), ExcDimensionMismatch(n(), dst.size()));
  Assert(m() == src.size(), ExcDimensionMismatch(m(), src.size()));

  Assert(!PointerComparison::equal(&src, &dst), ExcSourceEqualsDestination());

  dst = 0;

  for (size_type i = 0; i < m(); ++i)
    {
      for (size_type j = cols->rowstart[i]; j < cols->rowstart[i + 1]; ++j)
        {
          const size_type p = cols->colnums[j];
          dst(p) += typename OutVector::value_type(val[j]) *
                    typename OutVector::value_type(src(i));
        }
    }
}



template <typename number>
template <class OutVector, class InVector>
void
SparseMatrix<number>::vmult_add(OutVector &dst, const InVector &src) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  Assert(m() == dst.size(), ExcDimensionMismatch(m(), dst.size()));
  Assert(n() == src.size(), ExcDimensionMismatch(n(), src.size()));

  Assert(!PointerComparison::equal(&src, &dst), ExcSourceEqualsDestination());

  parallel::apply_to_subranges(
    0U,
    m(),
    [this, &src, &dst](const size_type begin_row, const size_type end_row) {
      internal::SparseMatrixImplementation::vmult_on_subrange(
        begin_row,
        end_row,
        val.get(),
        cols->rowstart.get(),
        cols->colnums.get(),
        src,
        dst,
        true);
    },
    internal::SparseMatrixImplementation::minimum_parallel_grain_size);
}



template <typename number>
template <class OutVector, class InVector>
void
SparseMatrix<number>::Tvmult_add(OutVector &dst, const InVector &src) const
{
  Assert(val != nullptr, ExcNotInitialized());
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(n() == dst.size(), ExcDimensionMismatch(n(), dst.size()));
  Assert(m() == src.size(), ExcDimensionMismatch(m(), src.size()));

  Assert(!PointerComparison::equal(&src, &dst), ExcSourceEqualsDestination());

  for (size_type i = 0; i < m(); ++i)
    for (size_type j = cols->rowstart[i]; j < cols->rowstart[i + 1]; ++j)
      {
        const size_type p = cols->colnums[j];
        dst(p) += typename OutVector::value_type(val[j]) *
                  typename OutVector::value_type(src(i));
      }
}


namespace internal
{
  namespace SparseMatrixImplementation
  {
    /**
     * Perform a vmult using the SparseMatrix data structures, but only using
     * a subinterval for the row indices.
     *
     * In the sequential case, this function is called on all rows, in the
     * parallel case it may be called on a subrange, at the discretion of the
     * task scheduler.
     */
    template <typename number, typename InVector>
    typename InVector::value_type
    matrix_norm_sqr_on_subrange(const size_type    begin_row,
                                const size_type    end_row,
                                const number      *values,
                                const std::size_t *rowstart,
                                const size_type   *colnums,
                                const InVector    &v)
    {
      typename InVector::value_type norm_sqr = 0.;

      for (size_type i = begin_row; i < end_row; ++i)
        {
          typename InVector::value_type s = 0;
          for (size_type j = rowstart[i]; j < rowstart[i + 1]; ++j)
            s += typename InVector::value_type(values[j]) * v(colnums[j]);
          norm_sqr +=
            v(i) *
            numbers::NumberTraits<typename InVector::value_type>::conjugate(s);
        }
      return norm_sqr;
    }
  } // namespace SparseMatrixImplementation
} // namespace internal



template <typename number>
template <typename somenumber>
somenumber
SparseMatrix<number>::matrix_norm_square(const Vector<somenumber> &v) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  Assert(m() == v.size(), ExcDimensionMismatch(m(), v.size()));
  Assert(n() == v.size(), ExcDimensionMismatch(n(), v.size()));

  return parallel::accumulate_from_subranges<somenumber>(
    [this, &v](const size_type begin_row, const size_type end_row) {
      return internal::SparseMatrixImplementation::matrix_norm_sqr_on_subrange(
        begin_row,
        end_row,
        val.get(),
        cols->rowstart.get(),
        cols->colnums.get(),
        v);
    },
    0,
    m(),
    internal::SparseMatrixImplementation::minimum_parallel_grain_size);
}



namespace internal
{
  namespace SparseMatrixImplementation
  {
    /**
     * Perform a vmult using the SparseMatrix data structures, but only using
     * a subinterval for the row indices.
     *
     * In the sequential case, this function is called on all rows, in the
     * parallel case it may be called on a subrange, at the discretion of the
     * task scheduler.
     */
    template <typename number, typename InVector>
    typename InVector::value_type
    matrix_scalar_product_on_subrange(const size_type    begin_row,
                                      const size_type    end_row,
                                      const number      *values,
                                      const std::size_t *rowstart,
                                      const size_type   *colnums,
                                      const InVector    &u,
                                      const InVector    &v)
    {
      typename InVector::value_type norm_sqr = 0.;

      for (size_type i = begin_row; i < end_row; ++i)
        {
          typename InVector::value_type s = 0;
          for (size_type j = rowstart[i]; j < rowstart[i + 1]; ++j)
            s += typename InVector::value_type(values[j]) * v(colnums[j]);
          norm_sqr +=
            u(i) *
            numbers::NumberTraits<typename InVector::value_type>::conjugate(s);
        }
      return norm_sqr;
    }
  } // namespace SparseMatrixImplementation
} // namespace internal



template <typename number>
template <typename somenumber>
somenumber
SparseMatrix<number>::matrix_scalar_product(const Vector<somenumber> &u,
                                            const Vector<somenumber> &v) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  Assert(m() == u.size(), ExcDimensionMismatch(m(), u.size()));
  Assert(n() == v.size(), ExcDimensionMismatch(n(), v.size()));

  return parallel::accumulate_from_subranges<somenumber>(
    [this, &u, &v](const size_type begin_row, const size_type end_row) {
      return internal::SparseMatrixImplementation::
        matrix_scalar_product_on_subrange(begin_row,
                                          end_row,
                                          val.get(),
                                          cols->rowstart.get(),
                                          cols->colnums.get(),
                                          u,
                                          v);
    },
    0,
    m(),
    internal::SparseMatrixImplementation::minimum_parallel_grain_size);
}



template <typename number>
template <typename numberB, typename numberC>
void
SparseMatrix<number>::mmult(SparseMatrix<numberC>       &C,
                            const SparseMatrix<numberB> &B,
                            const Vector<number>        &V,
                            const bool rebuild_sparsity_C) const
{
  const bool use_vector = V.size() == n() ? true : false;
  Assert(n() == B.m(), ExcDimensionMismatch(n(), B.m()));
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(B.cols != nullptr, ExcNeedsSparsityPattern());
  Assert(C.cols != nullptr, ExcNeedsSparsityPattern());

  const SparsityPattern &sp_A = *cols;
  const SparsityPattern &sp_B = *B.cols;

  // clear previous content of C
  if (rebuild_sparsity_C == true)
    {
      // we are about to change the sparsity pattern of C. this can not work
      // if either A or B use the same sparsity pattern
      Assert(&C.get_sparsity_pattern() != &this->get_sparsity_pattern(),
             ExcMessage("Can't use the same sparsity pattern for "
                        "different matrices if it is to be rebuilt."));
      Assert(&C.get_sparsity_pattern() != &B.get_sparsity_pattern(),
             ExcMessage("Can't use the same sparsity pattern for "
                        "different matrices if it is to be rebuilt."));

      // need to change the sparsity pattern of C, so cast away const-ness.
      SparsityPattern &sp_C =
        *(const_cast<SparsityPattern *>(&C.get_sparsity_pattern()));
      C.clear();
      sp_C.reinit(0, 0, 0);

      // create a sparsity pattern for the matrix C.
      {
        DynamicSparsityPattern dsp;
        dsp.compute_mmult_pattern(sp_A, sp_B);
        sp_C.copy_from(dsp);
      }

      // reinit matrix C from that information
      C.reinit(sp_C);
    }

  Assert(C.m() == m(), ExcDimensionMismatch(C.m(), m()));
  Assert(C.n() == B.n(), ExcDimensionMismatch(C.n(), B.n()));

  // create an array that caches some
  // elements that are going to be written
  // into the new matrix.
  unsigned int max_n_cols_B = 0;
  for (size_type i = 0; i < B.m(); ++i)
    max_n_cols_B = std::max(max_n_cols_B, sp_B.row_length(i));
  std::vector<numberC> new_entries(max_n_cols_B);

  // now compute the actual entries: a matrix-matrix product involves three
  // nested loops. One over the rows of A, for each row we then loop over all
  // the columns, and then we need to multiply each element with all the
  // elements in that row in B.
  for (size_type i = 0; i < C.m(); ++i)
    {
      const size_type       *rows     = &sp_A.colnums[sp_A.rowstart[i]];
      const size_type *const end_rows = &sp_A.colnums[sp_A.rowstart[i + 1]];
      for (; rows != end_rows; ++rows)
        {
          const number     A_val = val[rows - &sp_A.colnums[sp_A.rowstart[0]]];
          const size_type  col   = *rows;
          const size_type *new_cols = (&sp_B.colnums[sp_B.rowstart[col]]);

          // special treatment for diagonal
          if (sp_B.n_rows() == sp_B.n_cols())
            {
              C.add(i,
                    *new_cols,
                    numberC(A_val) *
                      numberC(
                        B.val[new_cols - &sp_B.colnums[sp_B.rowstart[0]]]) *
                      numberC(use_vector ? V(col) : 1));
              ++new_cols;
            }

          // now the innermost loop that goes over all the elements in row
          // 'col' of matrix B. Cache the elements, and then write them into C
          // at once
          numberC       *new_ptr = new_entries.data();
          const numberB *B_val_ptr =
            &B.val[new_cols - &sp_B.colnums[sp_B.rowstart[0]]];
          const numberB *const end_cols = &B.val[sp_B.rowstart[col + 1]];
          for (; B_val_ptr != end_cols; ++B_val_ptr)
            *new_ptr++ = numberC(A_val) * numberC(*B_val_ptr) *
                         numberC(use_vector ? V(col) : 1);

          C.add(i,
                new_ptr - new_entries.data(),
                new_cols,
                new_entries.data(),
                false,
                true);
        }
    }
}



template <typename number>
template <typename numberB, typename numberC>
void
SparseMatrix<number>::Tmmult(SparseMatrix<numberC>       &C,
                             const SparseMatrix<numberB> &B,
                             const Vector<number>        &V,
                             const bool rebuild_sparsity_C) const
{
  const bool use_vector = V.size() == m() ? true : false;
  Assert(m() == B.m(), ExcDimensionMismatch(m(), B.m()));
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(B.cols != nullptr, ExcNeedsSparsityPattern());
  Assert(C.cols != nullptr, ExcNeedsSparsityPattern());

  const SparsityPattern &sp_A = *cols;
  const SparsityPattern &sp_B = *B.cols;

  // clear previous content of C
  if (rebuild_sparsity_C == true)
    {
      // we are about to change the sparsity pattern of C. this can not work
      // if either A or B use the same sparsity pattern
      Assert(&C.get_sparsity_pattern() != &this->get_sparsity_pattern(),
             ExcMessage("Can't use the same sparsity pattern for "
                        "different matrices if it is to be rebuilt."));
      Assert(&C.get_sparsity_pattern() != &B.get_sparsity_pattern(),
             ExcMessage("Can't use the same sparsity pattern for "
                        "different matrices if it is to be rebuilt."));

      // need to change the sparsity pattern of C, so cast away const-ness.
      SparsityPattern &sp_C =
        *(const_cast<SparsityPattern *>(&C.get_sparsity_pattern()));
      C.clear();
      sp_C.reinit(0, 0, 0);

      // create a sparsity pattern for the matrix.
      {
        DynamicSparsityPattern dsp;
        dsp.compute_Tmmult_pattern(sp_A, sp_B);
        sp_C.copy_from(dsp);
      }

      // reinit matrix C from that information
      C.reinit(sp_C);
    }

  Assert(C.m() == n(), ExcDimensionMismatch(C.m(), n()));
  Assert(C.n() == B.n(), ExcDimensionMismatch(C.n(), B.n()));

  // create an array that caches some
  // elements that are going to be written
  // into the new matrix.
  unsigned int max_n_cols_B = 0;
  for (size_type i = 0; i < B.m(); ++i)
    max_n_cols_B = std::max(max_n_cols_B, sp_B.row_length(i));
  std::vector<numberC> new_entries(max_n_cols_B);

  // now compute the actual entries: a matrix-matrix product involves three
  // nested loops. One over the rows of A, for each row we then loop over all
  // the columns, and then we need to multiply each element with all the
  // elements in that row in B.
  for (size_type i = 0; i < m(); ++i)
    {
      const size_type       *rows     = &sp_A.colnums[sp_A.rowstart[i]];
      const size_type *const end_rows = &sp_A.colnums[sp_A.rowstart[i + 1]];
      const size_type       *new_cols = &sp_B.colnums[sp_B.rowstart[i]];
      if (sp_B.n_rows() == sp_B.n_cols())
        ++new_cols;

      const numberB *const end_cols = &B.val[sp_B.rowstart[i + 1]];

      for (; rows != end_rows; ++rows)
        {
          const size_type row   = *rows;
          const number    A_val = val[rows - &sp_A.colnums[sp_A.rowstart[0]]];

          // special treatment for diagonal
          if (sp_B.n_rows() == sp_B.n_cols())
            C.add(row,
                  i,
                  numberC(A_val) *
                    numberC(
                      B.val[new_cols - 1 - &sp_B.colnums[sp_B.rowstart[0]]]) *
                    numberC(use_vector ? V(i) : 1));

          // now the innermost loop that goes over all the elements in row
          // 'col' of matrix B. Cache the elements, and then write them into C
          // at once
          numberC       *new_ptr = new_entries.data();
          const numberB *B_val_ptr =
            &B.val[new_cols - &sp_B.colnums[sp_B.rowstart[0]]];
          for (; B_val_ptr != end_cols; ++B_val_ptr)
            *new_ptr++ = numberC(A_val) * numberC(*B_val_ptr) *
                         numberC(use_vector ? V(i) : 1);

          C.add(row,
                new_ptr - new_entries.data(),
                new_cols,
                new_entries.data(),
                false,
                true);
        }
    }
}



template <typename number>
typename SparseMatrix<number>::real_type
SparseMatrix<number>::l1_norm() const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());

  Vector<real_type> column_sums(n());
  const size_type   n_rows = m();
  for (size_type row = 0; row < n_rows; ++row)
    for (size_type j = cols->rowstart[row]; j < cols->rowstart[row + 1]; ++j)
      column_sums(cols->colnums[j]) +=
        numbers::NumberTraits<number>::abs(val[j]);

  return column_sums.linfty_norm();
}



template <typename number>
typename SparseMatrix<number>::real_type
SparseMatrix<number>::linfty_norm() const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());

  const number *val_ptr = &val[cols->rowstart[0]];

  real_type       max    = 0;
  const size_type n_rows = m();
  for (size_type row = 0; row < n_rows; ++row)
    {
      real_type           sum            = 0;
      const number *const val_end_of_row = &val[cols->rowstart[row + 1]];
      while (val_ptr != val_end_of_row)
        sum += numbers::NumberTraits<number>::abs(*val_ptr++);
      if (sum > max)
        max = sum;
    }
  return max;
}



template <typename number>
typename SparseMatrix<number>::real_type
SparseMatrix<number>::frobenius_norm() const
{
  // simply add up all entries in the
  // sparsity pattern, without taking any
  // reference to rows or columns
  real_type       norm_sqr = 0;
  const size_type n_rows   = m();
  for (const number *ptr = val.get(); ptr != val.get() + cols->rowstart[n_rows];
       ++ptr)
    norm_sqr += numbers::NumberTraits<number>::abs_square(*ptr);

  return std::sqrt(norm_sqr);
}



namespace internal
{
  namespace SparseMatrixImplementation
  {
    /**
     * Perform a vmult using the SparseMatrix data structures, but only using
     * a subinterval for the row indices.
     *
     * In the sequential case, this function is called on all rows, in the
     * parallel case it may be called on a subrange, at the discretion of the
     * task scheduler.
     */
    template <typename number, typename InVector, typename OutVector>
    typename OutVector::value_type
    residual_sqr_on_subrange(const size_type    begin_row,
                             const size_type    end_row,
                             const number      *values,
                             const std::size_t *rowstart,
                             const size_type   *colnums,
                             const InVector    &u,
                             const InVector    &b,
                             OutVector         &dst)
    {
      typename OutVector::value_type norm_sqr = 0.;

      for (size_type i = begin_row; i < end_row; ++i)
        {
          typename OutVector::value_type s = b(i);
          for (size_type j = rowstart[i]; j < rowstart[i + 1]; ++j)
            s -= typename OutVector::value_type(values[j]) * u(colnums[j]);
          dst(i) = s;
          norm_sqr +=
            s *
            numbers::NumberTraits<typename OutVector::value_type>::conjugate(s);
        }
      return norm_sqr;
    }
  } // namespace SparseMatrixImplementation
} // namespace internal


template <typename number>
template <typename somenumber>
somenumber
SparseMatrix<number>::residual(Vector<somenumber>       &dst,
                               const Vector<somenumber> &u,
                               const Vector<somenumber> &b) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  Assert(m() == dst.size(), ExcDimensionMismatch(m(), dst.size()));
  Assert(m() == b.size(), ExcDimensionMismatch(m(), b.size()));
  Assert(n() == u.size(), ExcDimensionMismatch(n(), u.size()));

  Assert(&u != &dst, ExcSourceEqualsDestination());

  return std::sqrt(parallel::accumulate_from_subranges<somenumber>(
    [this, &u, &b, &dst](const size_type begin_row, const size_type end_row) {
      return internal::SparseMatrixImplementation::residual_sqr_on_subrange(
        begin_row,
        end_row,
        val.get(),
        cols->rowstart.get(),
        cols->colnums.get(),
        u,
        b,
        dst);
    },
    0,
    m(),
    internal::SparseMatrixImplementation::minimum_parallel_grain_size));
}


namespace internal
{
  namespace SparseMatrixImplementation
  {
    // assert that the matrix has no zeros on the diagonal. this is important
    // for preconditioners such as Jacobi or SSOR
    template <typename number>
    void
    AssertNoZerosOnDiagonal(const SparseMatrix<number> &matrix)
    {
      if constexpr (running_in_debug_mode())
        {
          for (typename SparseMatrix<number>::size_type row = 0;
               row < matrix.m();
               ++row)
            Assert(matrix.diag_element(row) != number(),
                   ExcMessage(
                     "There is a zero on the diagonal of this matrix "
                     "in row " +
                     std::to_string(row) +
                     ". The preconditioner you selected cannot work if that "
                     "is the case because one of its steps requires "
                     "division by the diagonal elements of the matrix."
                     "\n\n"
                     "You should check whether you have correctly "
                     "assembled the matrix that you use for this "
                     "preconditioner. If it is correct that there are "
                     "zeros on the diagonal, then you will have to chose "
                     "a different preconditioner."));
        }
      else
        {
          (void)matrix;
        }
    }
  } // namespace SparseMatrixImplementation
} // namespace internal


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::precondition_Jacobi(Vector<somenumber>       &dst,
                                          const Vector<somenumber> &src,
                                          const number              omega) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  AssertDimension(m(), n());
  AssertDimension(dst.size(), n());
  AssertDimension(src.size(), n());

  internal::SparseMatrixImplementation::AssertNoZerosOnDiagonal(*this);

  const size_type    n            = src.size();
  somenumber        *dst_ptr      = dst.begin();
  const somenumber  *src_ptr      = src.begin();
  const std::size_t *rowstart_ptr = cols->rowstart.get();

  // optimize the following loop for
  // the case that the relaxation
  // factor is one. In that case, we
  // can save one FP multiplication
  // per row
  //
  // note that for square matrices,
  // the diagonal entry is the first
  // in each row, i.e. at index
  // rowstart[i]. and we do have a
  // square matrix by above assertion
  if (omega != number(1.))
    for (size_type i = 0; i < n; ++i, ++dst_ptr, ++src_ptr, ++rowstart_ptr)
      *dst_ptr = somenumber(omega) * *src_ptr / somenumber(val[*rowstart_ptr]);
  else
    for (size_type i = 0; i < n; ++i, ++dst_ptr, ++src_ptr, ++rowstart_ptr)
      *dst_ptr = *src_ptr / somenumber(val[*rowstart_ptr]);
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::precondition_SSOR(
  Vector<somenumber>             &dst,
  const Vector<somenumber>       &src,
  const number                    omega,
  const std::vector<std::size_t> &pos_right_of_diagonal) const
{
  // to understand how this function works
  // you may want to take a look at the CVS
  // archives to see the original version
  // which is much clearer...
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  AssertDimension(m(), n());
  AssertDimension(dst.size(), n());
  AssertDimension(src.size(), n());

  internal::SparseMatrixImplementation::AssertNoZerosOnDiagonal(*this);

  const size_type    n            = src.size();
  const std::size_t *rowstart_ptr = cols->rowstart.get();
  somenumber        *dst_ptr      = &dst(0);

  // case when we have stored the position
  // just right of the diagonal (then we
  // don't have to search for it).
  if (pos_right_of_diagonal.size() != 0)
    {
      Assert(pos_right_of_diagonal.size() == dst.size(),
             ExcDimensionMismatch(pos_right_of_diagonal.size(), dst.size()));

      // forward sweep
      for (size_type row = 0; row < n; ++row, ++dst_ptr, ++rowstart_ptr)
        {
          *dst_ptr = src(row);
          const std::size_t first_right_of_diagonal_index =
            pos_right_of_diagonal[row];
          Assert(first_right_of_diagonal_index <= *(rowstart_ptr + 1),
                 ExcInternalError());
          number s = 0;
          for (size_type j = (*rowstart_ptr) + 1;
               j < first_right_of_diagonal_index;
               ++j)
            s += val[j] * number(dst(cols->colnums[j]));

          // divide by diagonal element
          *dst_ptr -= s * omega;
          *dst_ptr /= val[*rowstart_ptr];
        }

      rowstart_ptr = cols->rowstart.get();
      dst_ptr      = &dst(0);
      for (; rowstart_ptr != &cols->rowstart[n]; ++rowstart_ptr, ++dst_ptr)
        *dst_ptr *= somenumber(omega * (number(2.) - omega)) *
                    somenumber(val[*rowstart_ptr]);

      // backward sweep
      rowstart_ptr = &cols->rowstart[n - 1];
      dst_ptr      = &dst(n - 1);
      for (int row = n - 1; row >= 0; --row, --rowstart_ptr, --dst_ptr)
        {
          const size_type end_row = *(rowstart_ptr + 1);
          const size_type first_right_of_diagonal_index =
            pos_right_of_diagonal[row];
          number s = 0;
          // go through the column from the end towards the diagonal in order
          // to delay the use of the newly computed "dst" values on
          // out-of-order-execution hardware
          for (size_type j = end_row - 1; j >= first_right_of_diagonal_index;
               --j)
            s += val[j] * number(dst(cols->colnums[j]));

          *dst_ptr -= s * omega;
          *dst_ptr /= val[*rowstart_ptr];
        };
      return;
    }

  // case when we need to get the position
  // of the first element right of the
  // diagonal manually for each sweep.
  // forward sweep
  for (size_type row = 0; row < n; ++row, ++dst_ptr, ++rowstart_ptr)
    {
      *dst_ptr = src(row);
      // find the first element in this line
      // which is on the right of the diagonal.
      // we need to precondition with the
      // elements on the left only.
      // note: the first entry in each
      // line denotes the diagonal element,
      // which we need not check.
      const size_type first_right_of_diagonal_index =
        (Utilities::lower_bound(cols->colnums.get() + *rowstart_ptr + 1,
                                cols->colnums.get() + *(rowstart_ptr + 1),
                                row) -
         cols->colnums.get());

      number s = 0;
      for (size_type j = (*rowstart_ptr) + 1; j < first_right_of_diagonal_index;
           ++j)
        s += val[j] * number(dst(cols->colnums[j]));

      // divide by diagonal element
      *dst_ptr -= s * omega;
      Assert(val[*rowstart_ptr] != number(), ExcDivideByZero());
      *dst_ptr /= val[*rowstart_ptr];
    };

  rowstart_ptr = cols->rowstart.get();
  dst_ptr      = &dst(0);
  for (size_type row = 0; row < n; ++row, ++rowstart_ptr, ++dst_ptr)
    *dst_ptr *=
      somenumber((number(2.) - omega)) * somenumber(val[*rowstart_ptr]);

  // backward sweep
  rowstart_ptr = &cols->rowstart[n - 1];
  dst_ptr      = &dst(n - 1);
  for (int row = n - 1; row >= 0; --row, --rowstart_ptr, --dst_ptr)
    {
      const size_type end_row = *(rowstart_ptr + 1);
      const size_type first_right_of_diagonal_index =
        (Utilities::lower_bound(&cols->colnums[*rowstart_ptr + 1],
                                &cols->colnums[end_row],
                                static_cast<size_type>(row)) -
         cols->colnums.get());
      number s = 0;
      for (size_type j = first_right_of_diagonal_index; j < end_row; ++j)
        s += val[j] * number(dst(cols->colnums[j]));
      *dst_ptr -= s * omega;
      Assert(val[*rowstart_ptr] != number(), ExcDivideByZero());
      *dst_ptr /= val[*rowstart_ptr];
    };
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::precondition_SOR(Vector<somenumber>       &dst,
                                       const Vector<somenumber> &src,
                                       const number              omega) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());

  dst = src;
  SOR(dst, omega);
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::precondition_TSOR(Vector<somenumber>       &dst,
                                        const Vector<somenumber> &src,
                                        const number              omega) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());

  dst = src;
  TSOR(dst, omega);
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::SOR(Vector<somenumber> &dst, const number omega) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  AssertDimension(m(), n());
  AssertDimension(dst.size(), n());

  internal::SparseMatrixImplementation::AssertNoZerosOnDiagonal(*this);

  for (size_type row = 0; row < m(); ++row)
    {
      somenumber s = dst(row);
      for (size_type j = cols->rowstart[row]; j < cols->rowstart[row + 1]; ++j)
        {
          const size_type col = cols->colnums[j];
          if (col < row)
            s -= somenumber(val[j]) * dst(col);
        }

      dst(row) = s * somenumber(omega) / somenumber(val[cols->rowstart[row]]);
    }
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::TSOR(Vector<somenumber> &dst, const number omega) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  AssertDimension(m(), n());
  AssertDimension(dst.size(), n());

  internal::SparseMatrixImplementation::AssertNoZerosOnDiagonal(*this);

  size_type row = m() - 1;
  while (true)
    {
      somenumber s = dst(row);
      for (size_type j = cols->rowstart[row]; j < cols->rowstart[row + 1]; ++j)
        if (cols->colnums[j] > row)
          s -= somenumber(val[j]) * dst(cols->colnums[j]);

      dst(row) = s * somenumber(omega) / somenumber(val[cols->rowstart[row]]);

      if (row == 0)
        break;

      --row;
    }
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::PSOR(Vector<somenumber>           &dst,
                           const std::vector<size_type> &permutation,
                           const std::vector<size_type> &inverse_permutation,
                           const number                  omega) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  AssertDimension(m(), n());

  Assert(m() == dst.size(), ExcDimensionMismatch(m(), dst.size()));
  Assert(m() == permutation.size(),
         ExcDimensionMismatch(m(), permutation.size()));
  Assert(m() == inverse_permutation.size(),
         ExcDimensionMismatch(m(), inverse_permutation.size()));

  internal::SparseMatrixImplementation::AssertNoZerosOnDiagonal(*this);

  for (size_type urow = 0; urow < m(); ++urow)
    {
      const size_type row = permutation[urow];
      somenumber      s   = dst(row);

      for (size_type j = cols->rowstart[row]; j < cols->rowstart[row + 1]; ++j)
        {
          const size_type col = cols->colnums[j];
          if (inverse_permutation[col] < urow)
            {
              s -= somenumber(val[j]) * dst(col);
            }
        }

      dst(row) = s * somenumber(omega) / somenumber(val[cols->rowstart[row]]);
    }
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::TPSOR(Vector<somenumber>           &dst,
                            const std::vector<size_type> &permutation,
                            const std::vector<size_type> &inverse_permutation,
                            const number                  omega) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  AssertDimension(m(), n());

  Assert(m() == dst.size(), ExcDimensionMismatch(m(), dst.size()));
  Assert(m() == permutation.size(),
         ExcDimensionMismatch(m(), permutation.size()));
  Assert(m() == inverse_permutation.size(),
         ExcDimensionMismatch(m(), inverse_permutation.size()));

  internal::SparseMatrixImplementation::AssertNoZerosOnDiagonal(*this);

  for (size_type urow = m(); urow != 0;)
    {
      --urow;
      const size_type row = permutation[urow];
      somenumber      s   = dst(row);
      for (size_type j = cols->rowstart[row]; j < cols->rowstart[row + 1]; ++j)
        {
          const size_type col = cols->colnums[j];
          if (inverse_permutation[col] > urow)
            s -= somenumber(val[j]) * dst(col);
        }

      dst(row) = s * somenumber(omega) / somenumber(val[cols->rowstart[row]]);
    }
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::Jacobi_step(Vector<somenumber>       &v,
                                  const Vector<somenumber> &b,
                                  const number              omega) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  AssertDimension(m(), n());

  Assert(m() == v.size(), ExcDimensionMismatch(m(), v.size()));
  Assert(m() == b.size(), ExcDimensionMismatch(m(), b.size()));

  GrowingVectorMemory<Vector<somenumber>>            mem;
  typename VectorMemory<Vector<somenumber>>::Pointer w(mem);
  w->reinit(v);

  if (!v.all_zero())
    {
      vmult(*w, v);
      *w -= b;
    }
  else
    w->equ(-1., b);
  precondition_Jacobi(*w, *w, omega);
  v -= *w;
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::SOR_step(Vector<somenumber>       &v,
                               const Vector<somenumber> &b,
                               const number              omega) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  AssertDimension(m(), n());
  Assert(m() == v.size(), ExcDimensionMismatch(m(), v.size()));
  Assert(m() == b.size(), ExcDimensionMismatch(m(), b.size()));

  internal::SparseMatrixImplementation::AssertNoZerosOnDiagonal(*this);

  for (size_type row = 0; row < m(); ++row)
    {
      somenumber s = b(row);
      for (size_type j = cols->rowstart[row]; j < cols->rowstart[row + 1]; ++j)
        {
          s -= somenumber(val[j]) * v(cols->colnums[j]);
        }
      v(row) += s * somenumber(omega) / somenumber(val[cols->rowstart[row]]);
    }
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::TSOR_step(Vector<somenumber>       &v,
                                const Vector<somenumber> &b,
                                const number              omega) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  AssertDimension(m(), n());
  Assert(m() == v.size(), ExcDimensionMismatch(m(), v.size()));
  Assert(m() == b.size(), ExcDimensionMismatch(m(), b.size()));

  internal::SparseMatrixImplementation::AssertNoZerosOnDiagonal(*this);

  for (int row = m() - 1; row >= 0; --row)
    {
      somenumber s = b(row);
      for (size_type j = cols->rowstart[row]; j < cols->rowstart[row + 1]; ++j)
        {
          s -= somenumber(val[j]) * v(cols->colnums[j]);
        }
      v(row) += s * somenumber(omega) / somenumber(val[cols->rowstart[row]]);
    }
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::SSOR_step(Vector<somenumber>       &v,
                                const Vector<somenumber> &b,
                                const number              omega) const
{
  SOR_step(v, b, omega);
  TSOR_step(v, b, omega);
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::SSOR(Vector<somenumber> &dst, const number omega) const
{
  // TODO: Is this called anywhere? If so, multiplication with omega(2-omega)D
  // is missing
  DEAL_II_NOT_IMPLEMENTED();

  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());
  AssertDimension(m(), n());
  Assert(m() == dst.size(), ExcDimensionMismatch(m(), dst.size()));

  internal::SparseMatrixImplementation::AssertNoZerosOnDiagonal(*this);

  const size_type n = dst.size();
  size_type       j;
  somenumber      s;

  for (size_type i = 0; i < n; ++i)
    {
      s = 0.;
      for (j = cols->rowstart[i]; j < cols->rowstart[i + 1]; ++j)
        {
          const size_type p = cols->colnums[j];
          if (p != SparsityPattern::invalid_entry)
            {
              if (i > j)
                s += somenumber(val[j]) * dst(p);
            }
        }
      dst(i) -= s * somenumber(omega);
      dst(i) /= somenumber(val[cols->rowstart[i]]);
    }

  for (int i = n - 1; i >= 0;
       i--) // this time, i is signed, but always positive!
    {
      s = 0.;
      for (j = cols->rowstart[i]; j < cols->rowstart[i + 1]; ++j)
        {
          const size_type p = cols->colnums[j];
          if (p != SparsityPattern::invalid_entry)
            {
              if (static_cast<size_type>(i) < j)
                s += somenumber(val[j]) * dst(p);
            }
        }
      dst(i) -= s * somenumber(omega) / somenumber(val[cols->rowstart[i]]);
    }
}



template <typename number>
void
SparseMatrix<number>::print_formatted(std::ostream      &out,
                                      const unsigned int precision,
                                      const bool         scientific,
                                      const unsigned int width_,
                                      const char        *zero_string,
                                      const double       denominator,
                                      const char        *separator) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());

  unsigned int width = width_;

  std::ios::fmtflags old_flags     = out.flags();
  unsigned int       old_precision = out.precision(precision);

  if (scientific)
    {
      out.setf(std::ios::scientific, std::ios::floatfield);
      if (width == 0u)
        width = precision + 7;
    }
  else
    {
      out.setf(std::ios::fixed, std::ios::floatfield);
      if (width == 0u)
        width = precision + 2;
    }

  for (size_type i = 0; i < m(); ++i)
    {
      for (size_type j = 0; j < n(); ++j)
        if ((*cols)(i, j) != SparsityPattern::invalid_entry)
          out << std::setw(width)
              << val[cols->operator()(i, j)] * number(denominator) << separator;
        else
          out << std::setw(width) << zero_string << separator;
      out << std::endl;
    };
  AssertThrow(out.fail() == false, ExcIO());

  // reset output format
  out.precision(old_precision);
  out.flags(old_flags);
}



template <typename number>
void
SparseMatrix<number>::print_pattern(std::ostream &out,
                                    const double  threshold) const
{
  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());

  for (size_type i = 0; i < m(); ++i)
    {
      for (size_type j = 0; j < n(); ++j)
        if ((*cols)(i, j) == SparsityPattern::invalid_entry)
          out << '.';
        else if (std::abs(val[cols->operator()(i, j)]) > threshold)
          out << '*';
        else
          out << ':';
      out << std::endl;
    };
  AssertThrow(out.fail() == false, ExcIO());
}



template <typename number>
void
SparseMatrix<number>::print_as_numpy_arrays(std::ostream      &out,
                                            const unsigned int precision) const
{
  AssertThrow(out.fail() == false, ExcIO());
  boost::io::ios_flags_saver restore_flags(out);

  out.precision(precision);

  Assert(cols != nullptr, ExcNeedsSparsityPattern());
  Assert(val != nullptr, ExcNotInitialized());

  std::vector<number> rows;
  std::vector<number> columns;
  std::vector<number> values;
  rows.reserve(n_nonzero_elements());
  columns.reserve(n_nonzero_elements());
  values.reserve(n_nonzero_elements());

  for (size_type i = 0; i < cols->rows; ++i)
    {
      for (size_type j = cols->rowstart[i]; j < cols->rowstart[i + 1]; ++j)
        {
          rows.push_back(i);
          columns.push_back(cols->colnums[j]);
          values.push_back(val[j]);
        }
    }

  for (auto d : values)
    out << d << ' ';
  out << '\n';

  for (auto r : rows)
    out << r << ' ';
  out << '\n';

  for (auto c : columns)
    out << c << ' ';
  out << '\n';
  out << std::flush;

  AssertThrow(out.fail() == false, ExcIO());
}



template <typename number>
void
SparseMatrix<number>::block_write(std::ostream &out) const
{
  AssertThrow(out.fail() == false, ExcIO());

  // first the simple objects,
  // bracketed in [...]
  out << '[' << max_len << "][";
  // then write out real data
  out.write(reinterpret_cast<const char *>(val.get()),
            reinterpret_cast<const char *>(val.get() + max_len) -
              reinterpret_cast<const char *>(val.get()));
  out << ']';

  AssertThrow(out.fail() == false, ExcIO());
}



template <typename number>
void
SparseMatrix<number>::block_read(std::istream &in)
{
  AssertThrow(in.fail() == false, ExcIO());

  char c;

  // first read in simple data
  in >> c;
  AssertThrow(c == '[', ExcIO());
  in >> max_len;

  in >> c;
  AssertThrow(c == ']', ExcIO());
  in >> c;
  AssertThrow(c == '[', ExcIO());

  // reallocate space
  val = std::make_unique<number[]>(max_len);

  // then read data
  in.read(reinterpret_cast<char *>(val.get()),
          reinterpret_cast<char *>(val.get() + max_len) -
            reinterpret_cast<char *>(val.get()));
  in >> c;
  AssertThrow(c == ']', ExcIO());
}


template <typename number>
void
SparseMatrix<number>::compress(VectorOperation::values)
{}


template <typename number>
std::size_t
SparseMatrix<number>::memory_consumption() const
{
  return max_len * static_cast<std::size_t>(sizeof(number)) + sizeof(*this);
}


DEAL_II_NAMESPACE_CLOSE

#endif
