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

#ifndef dealii_lapack_banded_matrix_h
#define dealii_lapack_banded_matrix_h

#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/lapack_support.h>

#include <memory>
#include <type_traits>

DEAL_II_NAMESPACE_OPEN

/**
 * Class encapsulating scalars stored by LAPACKBandedMatrix.
 */
class LAPACKBandedMatrixData
{
public:
  /**
   * Size type.
   */
  typedef std::make_unsigned<types::blas_int>::type size_type;

  /**
   * Default constructor.
   */
  LAPACKBandedMatrixData();

  /**
   * Constructor.
   */
  LAPACKBandedMatrixData(const size_type n_rows,
                         const size_type n_cols,
                         const size_type n_subdiagonals,
                         const size_type n_superdiagonals);

protected:
  /**
   * Number of rows in the matrix.
   */
  types::blas_int n_rows;

  /**
   * Number of columns in the matrix.
   */
  types::blas_int n_cols;

  /**
   * Number of subdiagonals (that is, diagonals below the main diagonal).
   */
  types::blas_int n_subdiagonals;

  /**
   * Number of superdiagonals (that is, diagonals above the main diagonal).
   */
  types::blas_int n_superdiagonals;

  /**
   * Number of superdiagonals actually stored by the matrix. This is larger
   * than the number of subdiagonals since the LU factorization routines
   * require additional space in the matrix.
   */
  types::blas_int n_allocated_superdiagonals;

  /**
   * LAPACKBandedMatrix stores, in a sense, two different matrices: the
   * logically banded matrix and the underlying two-dimensional array where
   * each row is a diagonal of the banded matrix. This parameter stores the
   * number of rows that are allocated in the second representation. Both
   * matrices have the same number of columns.
   */
  types::blas_int n_allocated_rows;
};

/**
 * A variant of LAPACKFullMatrix that uses a different storage scheme
 * optimized for matrices with only a small number of subdiagonals and
 * superdiagonals (i.e., a banded matrix). For example, here is a $5\times 5$
 * band matrix with two subdiagonals and two superdiagonals:
 *
 * \f[
 * \begin{pmatrix}
 * a00 & a01                   \\
 * a10 & a11 & a12             \\
 * a20 & a21 & a22 & a23       \\
 * {}  & a31 & a32 & a33 & a34 \\
 * {}  & {}  & a42 & a43 & a44 \\
 * \end{pmatrix}
 * \f]
 *
 * LAPACK defines a special storage scheme for banded matrices: essentially,
 * each diagonal of the banded matrix is stored as one row in a rectangular
 * matrix. This class stores the banded matrix in the correct format for
 * LAPACK but has essentially the same dense matrix interface as FullMatrix
 * and LAPACKFullMatrix.
 *
 * The primary advantage of this class is that it provides very accurate and
 * efficient linear solvers that take advantage of the band structure. It is
 * conceptually similar to SparseMatrix in the sense that it only stores the
 * nonzero entries (as well as a small amount of padding) of the matrix it
 * represents.
 *
 * @ingroup Matrix1
 * @author David Wells, 2018
 */
template <typename Number>
class LAPACKBandedMatrix : public Subscriptor, protected LAPACKBandedMatrixData
{
public:
  /**
   * Declare type for container size.
   */
  typedef std::make_unsigned<types::blas_int>::type size_type;

  /**
   * Default constructor. Creates a size zero empty matrix.
   */
  LAPACKBandedMatrix() = default;

  /**
   * Constructor. Creates a logically @p n_rows by @p n_cols matrix that
   * stores @p n_subdiagonals subdiagonals and @p n_superdiagonals
   * superdiagonals. The buffer stored by this array is in the format specified
   * by LAPACK: see
   *
   * http://www.netlib.org/lapack/lug/node124.html
   *
   * for a complete description of the format.
   */
  LAPACKBandedMatrix(const size_type n_rows,
                     const size_type n_cols,
                     const size_type n_subdiagonals,
                     const size_type n_superdiagonals);

  /**
   * Reinitialize the matrix to have the same state as if it were created by
   * the constructor of this class with the same arguments as the present
   * function.
   */
  void
  reinit(const size_type n_rows,
         const size_type n_cols,
         const size_type n_subdiagonals,
         const size_type n_superdiagonals);

  /**
   * Move constructor.
   */
  LAPACKBandedMatrix(LAPACKBandedMatrix<Number> &&) = default;

  /**
   * Move assignment operator.
   */
  LAPACKBandedMatrix<Number> &
  operator=(LAPACKBandedMatrix<Number> &&) = default;

  /**
   * Exchange the stored values of the current object with those in @p other.
   */
  void
  swap(LAPACKBandedMatrix<Number> &other);

  /**
   * Return the number of subdiagonals stored by the matrix.
   */
  size_type
  n_stored_subdiagonals() const;

  /**
   * Return the number of superdiagonals stored by the matrix.
   */
  size_type
  n_stored_superdiagonals() const;

  /**
   * Return the dimension of the codomain (or range) space.
   *
   * @note The matrix is of dimension $m \times n$.
   */
  size_type
  m() const;

  /**
   * Return the dimension of the domain space.
   *
   * @note The matrix is of dimension $m \times n$.
   */
  size_type
  n() const;

  /**
   * Return whether or not a given entry is stored in the matrix band. This
   * function is analogous to SparsityPattern::exists() and is only valid when
   * <code>row_n < m()</code> and <code>col_n < n()</code>.
   */
  bool
  exists(const size_type row_n, const size_type col_n) const;

  /**
   * Constant access operator for an entry in one of the matrix bands.
   */
  const Number &
  operator()(const size_type row_n, const size_type col_n) const;

  /**
   * Access operator for an entry in one of the matrix bands.
   */
  Number &
  operator()(const size_type row_n, const size_type col_n);

  /**
   * General access operator. If the requested entry is outside of the stored
   * matrix bands then this function returns zero.
   */
  Number
  el(const size_type row_n, const size_type col_n) const;

  /**
   * Set the element (@p row_n, @p col_n) to <tt>value</tt>.
   *
   * Like el, this function will throw an error if both the entry is not
   * stored and @p value is nonzero. It will also throw an error if @p value
   * is not finite.
   */
  void
  set(const size_type row_n, const size_type col_n, const Number value);

protected:
  /**
   * Array containing the actual values of the matrix in column-major
   * order. This array corresponds to a logically <code>(2*n_subdiagonals + 1
   * + n_subdiagonals)</code> by <code>n_cols</code> dense matrix. The first
   * <code>n_subdiagonals*n_cols<code> are only used during LU factorization
   * and do not correspond to values of the banded matrix.
   */
  AlignedVector<Number> values;

  /**
   * Array of pivot indices used by compute_lu_factorization and solve.
   */
  std::vector<types::blas_int> pivots;

  /**
   * Enumeration describing the current state of the matrix; e.g., if
   * compute_lu_factorization() is called then the original matrix is replaced
   * by the LU factorization and the state is changed from 'matrix' to 'lu'.
   */
  LAPACKSupport::State state;

  /**
   * Thread mutex.
   */
  mutable Threads::Mutex mutex;

  /**
   * Numerical work array.
   */
  mutable std::vector<Number> work;

  /**
   * Index work array.
   */
  mutable std::vector<types::blas_int> index_work;

  /**
   * LAPACK routines will overwrite the right-hand side with the solution:
   * however, the solution postprocessor requires the original right-hand
   * side. Store a temporary copy here.
   */
  mutable Vector<Number> temporary_solution;

  /**
   * LU factorization overwrites the content of the matrix: however, the
   * solution postprocessor requires the original matrix. Get around this by
   * storing a copy of the original matrix if we LU-factorize.
   */
  std::unique_ptr<LAPACKBandedMatrix<Number>> original_matrix;

private:
  /**
   * Auxiliary function that computes an index into LAPACKBandedMatrix::values
   * for a given entry.
   */
  std::size_t
  compute_index(const size_type row_n, const size_type col_n) const;
};



#ifndef DOXYGEN
template <typename Number>
inline const Number &
LAPACKBandedMatrix<Number>::operator()(const size_type row_n,
                                       const size_type col_n) const
{
  return values[compute_index(row_n, col_n)];
}



template <typename Number>
inline Number &
LAPACKBandedMatrix<Number>::operator()(const size_type row_n,
                                       const size_type col_n)
{
  return values[compute_index(row_n, col_n)];
}



template <typename Number>
inline Number
LAPACKBandedMatrix<Number>::el(const size_type row_n,
                               const size_type col_n) const
{
  if (exists(row_n, col_n))
    return operator()(row_n, col_n);
  return Number();
}



template <typename Number>
inline void
LAPACKBandedMatrix<Number>::set(const size_type row_n,
                                const size_type col_n,
                                const Number    value)
{
  AssertIsFinite(value);
  // permit setting nonexistent entries to zero
  if (value == Number() && !exists(row_n, col_n))
    return;
  operator()(row_n, col_n) = value;
}



template <typename Number>
inline typename LAPACKBandedMatrix<Number>::size_type
LAPACKBandedMatrix<Number>::n_stored_subdiagonals() const
{
  return static_cast<LAPACKBandedMatrixData::size_type>(n_subdiagonals);
}



template <typename Number>
inline typename LAPACKBandedMatrix<Number>::size_type
LAPACKBandedMatrix<Number>::n_stored_superdiagonals() const
{
  return static_cast<LAPACKBandedMatrixData::size_type>(n_superdiagonals);
}



template <typename Number>
inline typename LAPACKBandedMatrix<Number>::size_type
LAPACKBandedMatrix<Number>::m() const
{
  return n_rows;
}



template <typename Number>
inline typename LAPACKBandedMatrix<Number>::size_type
LAPACKBandedMatrix<Number>::n() const
{
  return n_cols;
}



template <typename Number>
inline bool
LAPACKBandedMatrix<Number>::exists(const size_type row_n,
                                   const size_type col_n) const
{
  Assert(row_n < m(), ExcIndexRange(row_n, 0, m()));
  Assert(col_n < n(), ExcIndexRange(col_n, 0, n()));

  if (row_n == col_n)
    return true;
  if (col_n < row_n) // below the main diagonal
    return static_cast<types::blas_int>(row_n - col_n) <= n_subdiagonals;
  if (row_n < col_n) // above the main diagonal
    return static_cast<types::blas_int>(col_n - row_n) <= n_superdiagonals;

  // we should not get here
  Assert(false, ExcInternalError());
  return false;
}



template <typename Number>
inline std::size_t
LAPACKBandedMatrix<Number>::compute_index(const size_type row_n,
                                          const size_type col_n) const
{
  Assert(exists(row_n, col_n),
         ExcMessage("The given entry (" + std::to_string(row_n) + ", " +
                    std::to_string(col_n) +
                    ") is not stored in one of the "
                    "matrix bands."));
  const size_type access_row_n = n_allocated_superdiagonals + row_n - col_n;
  const size_type access_col_n = col_n;
  return access_col_n * n_allocated_rows + access_row_n;
}
#endif

DEAL_II_NAMESPACE_CLOSE

#endif
