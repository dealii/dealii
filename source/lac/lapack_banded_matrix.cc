// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

#include <deal.II/base/std_cxx14/memory.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_banded_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/utilities.h>
#include <deal.II/lac/vector.h>

DEAL_II_NAMESPACE_OPEN

// ----------------------------------------------------------------------------
// LAPACKBandedMatrixData
// ----------------------------------------------------------------------------

LAPACKBandedMatrixData::LAPACKBandedMatrixData() :
  LAPACKBandedMatrixData(0, 0, 0, 0)
{}



LAPACKBandedMatrixData::LAPACKBandedMatrixData(
  const size_type n_rows,
  const size_type n_cols,
  const size_type n_subdiagonals,
  const size_type n_superdiagonals) :
  n_rows(n_rows),
  n_cols(n_cols),
  n_subdiagonals(n_subdiagonals),
  n_superdiagonals(n_superdiagonals),
  n_allocated_superdiagonals(n_superdiagonals + n_subdiagonals),
  n_allocated_rows(n_allocated_superdiagonals + n_subdiagonals + 1)
{}

// ----------------------------------------------------------------------------
// construction, copying, swapping, and reinitialization
// ----------------------------------------------------------------------------

template <typename Number>
LAPACKBandedMatrix<Number>::LAPACKBandedMatrix(
  const size_type n_rows,
  const size_type n_cols,
  const size_type n_subdiagonals,
  const size_type n_superdiagonals) :
  LAPACKBandedMatrixData(n_rows, n_cols, n_subdiagonals, n_superdiagonals)
{
  values.resize(n_allocated_rows * n_cols);
}



template <typename Number>
LAPACKBandedMatrix<Number>::LAPACKBandedMatrix(
  const LAPACKBandedMatrix<Number> &other) :
  Subscriptor(),
  LAPACKBandedMatrixData(other),
  values(other.values),
  pivots(other.pivots),
  state(other.state)
{
  if (other.original_matrix != nullptr)
    original_matrix = std_cxx14::make_unique<LAPACKBandedMatrix<Number>>(
      *other.original_matrix);
}



template <typename Number>
void
LAPACKBandedMatrix<Number>::reinit(const size_type n_rows,
                                   const size_type n_cols,
                                   const size_type n_subdiagonals,
                                   const size_type n_superdiagonals)
{
  LAPACKBandedMatrix<Number> other(
    n_rows, n_cols, n_subdiagonals, n_superdiagonals);
  swap(other);
}


namespace internal
{
  namespace LAPACKBandedMatrixImplementation
  {
    /*
     * Copy entries with the generic iterator interface.
     *
     * @warning This function assumes that the iterators will only point to
     * entries that are actually stored in @p matrix.
     */
    template <typename Iterator, typename Number>
    void
    iterator_copy(Iterator                    entry,
                  const Iterator              end,
                  const bool                  elide_zeros,
                  LAPACKBandedMatrix<Number> &matrix)
    {
      // TODO is there a good way to add a static_assert here, given that the
      // relevant iterators all return proxy objects?
      for (; entry != end; ++entry)
        {
          if (elide_zeros && entry->value() == Number(0.0))
            continue;
          matrix(entry->row(), entry->column()) = entry->value();
        }
    }

    /**
     * Copy a matrix (e.g., FullMatrix or LAPACKFullMatrix) into a banded
     * matrix.
     */
    template <typename MatrixType, typename Number>
    void
    matrix_copy(const MatrixType &          other,
                const bool                  elide_zeros,
                LAPACKBandedMatrix<Number> &matrix)
    {
      static_assert(
        std::is_same<typename MatrixType::value_type, Number>::value,
        "This copy iterator is only implemented for identical "
        "numeric types.");
      std::size_t n_subdiagonals   = 0;
      std::size_t n_superdiagonals = 0;

      for (const auto entry : other)
        {
          if (entry.value() == 0.0 && elide_zeros)
            continue;
          else
            {
              if (entry.column() < entry.row())
                // row index is higher so we are on a subdiagonal
                n_subdiagonals = std::max<std::size_t>(
                  n_subdiagonals, entry.row() - entry.column());
              else if (entry.row() < entry.column())
                // column index is higher so we are on a superdiagonal
                n_superdiagonals = std::max<std::size_t>(
                  n_superdiagonals, entry.column() - entry.row());
            }
        }
      matrix.reinit(other.m(), other.n(), n_subdiagonals, n_superdiagonals);
      internal::LAPACKBandedMatrixImplementation::iterator_copy(
        other.begin(), other.end(), elide_zeros, matrix);
    }
  } // namespace LAPACKBandedMatrixImplementation
} // namespace internal



template <typename Number>
LAPACKBandedMatrix<Number>::LAPACKBandedMatrix(
  const LAPACKFullMatrix<Number> &other)
{
  internal::LAPACKBandedMatrixImplementation::matrix_copy(other, true, *this);
}



template <typename Number>
LAPACKBandedMatrix<Number>::LAPACKBandedMatrix(const FullMatrix<Number> &other)
{
  internal::LAPACKBandedMatrixImplementation::matrix_copy(other, true, *this);
}



template <typename Number>
void
LAPACKBandedMatrix<Number>::swap(LAPACKBandedMatrix<Number> &other)
{
  std::swap(*this, other);
}



template <typename Number>
LAPACKBandedMatrix<Number>::LAPACKBandedMatrix(const SparsityPattern &other)
{
  // A SparsityPattern knows its bandwidth, but does not know the number of
  // subdiagonals and superdiagonals: calculate these here
  std::size_t n_other_subdiagonals   = 0;
  std::size_t n_other_superdiagonals = 0;

  for (const SparsityPatternIterators::Accessor entry : other)
    {
      if (entry.column() < entry.row())
        // row index is higher so we are on a subdiagonal
        n_other_subdiagonals = std::max<std::size_t>(
          n_other_subdiagonals, entry.row() - entry.column());
      else if (entry.row() < entry.column())
        // column index is higher so we are on a superdiagonal
        n_other_superdiagonals = std::max<std::size_t>(
          n_other_superdiagonals, entry.column() - entry.row());
    }
  this->reinit(other.n_rows(),
               other.n_cols(),
               n_other_subdiagonals,
               n_other_superdiagonals);
}



template <typename Number>
LAPACKBandedMatrix<Number>::LAPACKBandedMatrix(
  const SparseMatrix<Number> &other) :
  LAPACKBandedMatrix(other.get_sparsity_pattern())
{
  internal::LAPACKBandedMatrixImplementation::iterator_copy(
    other.begin(), other.end(), true, *this);
}



template <typename Number>
LAPACKBandedMatrix<Number>::LAPACKBandedMatrix(
  const SparseMatrixEZ<Number> &other)
{
  internal::LAPACKBandedMatrixImplementation::matrix_copy(other,
                                                          /*elide_zeros=*/false,
                                                          *this);
}


// ----------------------------------------------------------------------------
// matrix-vector multiplication
// ----------------------------------------------------------------------------

template <typename Number>
void
LAPACKBandedMatrix<Number>::vmult(Vector<Number> &      w,
                                  const Vector<Number> &v,
                                  const bool            adding) const
{
  do_vmult(w, v, LAPACKSupport::N, adding);
}



template <typename Number>
void
LAPACKBandedMatrix<Number>::vmult_add(Vector<Number> &      w,
                                      const Vector<Number> &v) const
{
  do_vmult(w, v, LAPACKSupport::N, true);
}



template <typename Number>
void
LAPACKBandedMatrix<Number>::Tvmult(Vector<Number> &      w,
                                   const Vector<Number> &v,
                                   const bool            adding) const
{
  do_vmult(w, v, LAPACKSupport::T, adding);
}



template <typename Number>
void
LAPACKBandedMatrix<Number>::Tvmult_add(Vector<Number> &      w,
                                       const Vector<Number> &v) const
{
  do_vmult(w, v, LAPACKSupport::T, true);
}



template <typename Number>
void
LAPACKBandedMatrix<Number>::do_vmult(Vector<Number> &      w,
                                     const Vector<Number> &v,
                                     const char            transpose,
                                     const bool            adding) const
{
  const Number alpha = 1.;
  const Number beta  = (adding ? 1. : 0.);

  // This is a bit tricky: we want to skip the first few rows in the matrix
  // (that correspond to diagonals only used by the LU factorization). To do
  // this we set LDA equal to the number of allocated rows and then skip the
  // first n_subdiagonals rows below when calling gbmv. This skips the correct
  // number of entries between each column.
  gbmv(&transpose,              // TRANS
       &n_rows,                 // M
       &n_cols,                 // N
       &n_subdiagonals,         // KL
       &n_superdiagonals,       // KU
       &alpha,                  // ALPHA
       &values[n_subdiagonals], // A
       &n_allocated_rows,       // LDA
       v.begin(),               // X
       &LAPACKSupport::one,     // INCX
       &beta,                   // BETA
       w.begin(),               // Y
       &LAPACKSupport::one      // INCY
  );
}



#include "lapack_banded_matrix.inst"


DEAL_II_NAMESPACE_CLOSE
