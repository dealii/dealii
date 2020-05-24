// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#ifndef dealii_qr_h
#define dealii_qr_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/utilities.h>

#include <boost/signals2/signal.hpp>

#include <memory>

DEAL_II_NAMESPACE_OPEN

/**
 * A base class for thin QR implementations.
 *
 * This class and classes derived from it are meant to build $Q$ and $R$
 * matrices one row/column at a time, i.e., by growing $R$ matrix from an empty
 * $0\times 0$ matrix to $N\times N$, where $N$ is the number of added column
 * vectors.
 *
 * As a consequence, matrices which have the same number of rows as each vector
 * (i.e. $Q$ matrix) is stored as a collection of vectors of `VectorType`.
 */
template <typename VectorType>
class BaseQR
{
  /**
   * Number type for R matrix.
   */
  using Number = typename VectorType::value_type;

protected:
  /**
   * Default private constructor.
   */
  BaseQR();

public:
  /**
   * Destructor.
   */
  virtual ~BaseQR() = default;

  /**
   * Append @p column to the QR factorization.
   * Returns <code>true</code> if the result is successful, i.e.
   * the columns are linearly independent. Otherwise the @p column
   * is rejected and the return value is <code>false</code>.
   */
  virtual bool
  append_column(const VectorType &column) = 0;

  /**
   * Remove a column @p k and update QR factorization.
   */
  virtual void
  remove_column(const unsigned int k = 0) = 0;

  /**
   * Return size of the subspace.
   */
  unsigned int
  size() const;

  /**
   * Return the current upper triangular matrix R.
   */
  const LAPACKFullMatrix<Number> &
  get_R() const;

  /**
   * Solve $Rx=y$. Vectors @p x and @p y should be consistent
   * with the current size of the subspace.
   * If @p transpose is <code>true</code>, $R^Tx=y$ is solved instead.
   */
  void
  solve(Vector<Number> &      x,
        const Vector<Number> &y,
        const bool            transpose = false) const;

  /**
   * Set $y = Qx$. The size of $x$ should be consistent with the
   * size of the R matrix.
   */
  virtual void
  multiply_with_Q(VectorType &y, const Vector<Number> &x) const = 0;

  /**
   * Set $y = Q^Tx$. The size of $x$ should be consistent with the
   * size of column vectors.
   */
  virtual void
  multiply_with_QT(Vector<Number> &y, const VectorType &x) const = 0;

  /**
   * Set $y = QRx$. The size of $x$ should be consistent with the
   * size of the R matrix.
   */
  virtual void
  multiply_with_A(VectorType &y, const Vector<Number> &x) const = 0;

  /**
   * Set $y = R^T Q^Tx$. The size of $x$ should be consistent with the
   * size of column vectors.
   */
  virtual void
  multiply_with_AT(Vector<Number> &y, const VectorType &x) const = 0;

  /**
   * Connect a slot to retrieve a notification when the Givens rotations
   * are performed.
   *
   * The function takes two indices, @p i and @p j, describing the plane of
   * rotation, and a triplet of numbers @p csr (cosine, sine and radius, see
   * Utilities::LinearAlgebra::givens_rotation()) which represents the rotation
   * matrix.
   */
  boost::signals2::connection
  connect_givens_slot(
    const std::function<void(const unsigned int           i,
                             const unsigned int           j,
                             const std::array<Number, 3> &csr)> &slot);

protected:
  /**
   * Compute $y=Hx$ where $H$ is the matrix formed by the column vectors stored
   * by this object.
   */
  void
  multiply_with_cols(VectorType &y, const Vector<Number> &x) const;

  /**
   * Multiply with transpose columns stored in the object.
   */
  void
  multiply_with_colsT(Vector<Number> &y, const VectorType &x) const;

  /**
   * A vector of unique pointers to store columns.
   */
  std::vector<std::unique_ptr<VectorType>> columns;

  /**
   * Matrix to store R.
   */
  LAPACKFullMatrix<Number> R;

  /**
   * Current size (number of columns in Q).
   */
  unsigned int current_size;

  /**
   * Signal used to retrieve a notification
   * when Givens rotations are performed in the `(i,j)`-plane.
   */
  boost::signals2::signal<void(const unsigned int i,
                               const unsigned int j,
                               const std::array<Number, 3> &)>
    givens_signal;
};

// clang-format off
/**
 * A class to compute and store the QR factorization of a matrix represented by a set of column vectors.
 *
 * The class is design to update a given (possibly empty) QR factorization
 * of a matrix $A$ (constructed incrementally by providing its columns)
 * due to the addition of a new column vector to $A$. This is equivalent to
 * constructing an orthonormal basis by the Gram-Schmidt procedure.
 * The class also provides update functionality when the first column
 * is removed.
 *
 * The `VectorType` template argument may either be a parallel and serial vector, and only need
 * to have basic operations such as additions, scalar product, etc.
 * It also needs to have a copy-constructor.
 *
 * See sections 6.5.2-6.5.3 on pp. 335-338 in
 * @code{.bib}
 * @Book{Golub2013,
 *   title     = {Matrix computations},
 *   publisher = {Johns Hopkins University Press},
 *   year      = {2013},
 *   author    = {Golub, Gene H and Van Loan, Charles F},
 *   edition   = {4},
 *  }
 * @endcode
 * as well as
 * @code{.bib}
 * @article{Daniel1976,
 *   author   = {Daniel, James W and Gragg, Walter Bill and Kaufman, Linda and Stewart, Gilbert W},
 *   title    = {{Reorthogonalization and stable algorithms for updating the Gram-Schmidt QR factorization}},
 *   journal  = {Mathematics of Computation},
 *   year     = {1976},
 *   volume   = {30},
 *   number   = {136},
 *   pages    = {772--795},
 * }
 * @Article{Reichel1990,
 *   author     = {Reichel, L. and Gragg, W. B.},
 *   title      = {{Algorithm 686: FORTRAN Subroutines for Updating the QR Decomposition}},
 *   journal    = {ACM Trans. Math. Softw.},
 *   year       = {1990},
 *   volume     = {16},
 *   number     = {4},
 *   pages      = {369--377},
 *   month      = dec,
 *   issn       = {0098-3500},
 *   acmid      = {98291},
 *   address    = {New York, NY, USA},
 *   doi        = {10.1145/98267.98291},
 *   issue_date = {Dec. 1990},
 *   numpages   = {9},
 *   publisher  = {ACM},
 *   url        = {http://doi.acm.org/10.1145/98267.98291},
 *  }
 * @endcode
 */
// clang-format on
template <typename VectorType>
class QR : public BaseQR<VectorType>
{
public:
  /**
   * Number type for R matrix.
   */
  using Number = typename VectorType::value_type;

  /**
   * Default constructor.
   */
  QR();

  /**
   * Destructor.
   */
  virtual ~QR() = default;

  /**
   * @copydoc BaseQR::append_column
   *
   * @note Currently this function always returns <code>true</code>.
   */
  virtual bool
  append_column(const VectorType &column);

  /**
   * Remove first column and update QR factorization.
   *
   * Starting from the given QR decomposition
   * $QR= A = [a_1\,\dots a_n], \quad a_i \in {\mathbb R}^m$
   * we aim at computing factorization of
   * $\tilde Q \tilde R= \tilde A = [a_2\,\dots a_n], \quad a_i \in {\mathbb
   * R}^m$.
   *
   * The standard approach is to partition $R$ as
   * \f[
   * R =
   * \begin{bmatrix}
   * r_{11} & w^T \\
   * 0      & R_{33}
   * \end{bmatrix}
   * \f]
   * It then follows that
   * \f[
   * Q^T \tilde A =
   * \begin{bmatrix}
   * 0 & w^T \\
   * 0 & R_{33}
   * \end{bmatrix}
   * \f]
   * is upper Hessenberg where unwanted sub-diagonal elements can be
   * zeroed by a sequence of Givens rotations.
   *
   * Note that $\tilde R^T \tilde R = \tilde A^T \tilde A$,
   * where the RHS is included in $A^T A = R^T R$. Therefore
   * $\tilde R$ can be obtained by Cholesky decomposition.
   */
  virtual void
  remove_column(const unsigned int k = 0);

  virtual void
  multiply_with_Q(VectorType &y, const Vector<Number> &x) const;

  virtual void
  multiply_with_QT(Vector<Number> &y, const VectorType &x) const;

  virtual void
  multiply_with_A(VectorType &y, const Vector<Number> &x) const;

  virtual void
  multiply_with_AT(Vector<Number> &y, const VectorType &x) const;

private:
  /**
   * Apply givens rotation in the `(i,j)`-plane to @p Q and @p R so that
   * <code>R(k,k)</code> is zeroed.
   *
   * See Chapter 5.1.9 of Golub 2013, Matrix computations.
   */
  void
  apply_givens_rotation(const unsigned int i, const unsigned int k);

  /**
   * Temporary vector needed to do Givens rotation of Q.
   */
  VectorType tmp;
};



/**
 * A class to obtain the triangular $R$ matrix of the $A=QR$ factorization
 * together with the matrix $A$ itself. The orthonormal matrix $Q$ is not stored
 * explicitly, the name of the class.
 * The multiplication with $Q$ can be represented as $Q=A R^{-1}$, whereas the
 * multiplication with $Q^T$ is given by $Q^T=R^{-T}A^T$.
 *
 * The class is designed to update a given (possibly empty) QR factorization
 * due to the addition of a new column vector. This is equivalent to
 * constructing an orthonormal basis by the Gram-Schmidt procedure.
 * The class also provides update functionality when the column
 * is removed.
 *
 * The `VectorType` template argument may either be a parallel and serial
 * vector, and only need to have basic operations such as additions, scalar
 * product, etc. It also needs to have a copy-constructor.
 */
template <typename VectorType>
class ImplicitQR : public BaseQR<VectorType>
{
public:
  /**
   * Number type for R matrix.
   */
  using Number = typename VectorType::value_type;

  /**
   * Default constructor.
   */
  ImplicitQR();

  /**
   * Destructor.
   */
  virtual ~ImplicitQR() = default;

  virtual bool
  append_column(const VectorType &column);

  /**
   * Remove column and update QR factorization.
   *
   * Starting from the given QR decomposition
   * $QR= A = [a_1\,\dots a_n], \quad a_i \in R^m$
   * we aim at computing factorization of
   * $\tilde Q \tilde R= \tilde A = [a_2\,\dots a_n], \quad a_i \in R^m$.
   *
   * Note that $\tilde R^T \tilde R = \tilde A^T \tilde A$,
   * where the RHS is included in $A^T A = R^T R$. Therefore
   * $\tilde R$ can be obtained by Cholesky decomposition.
   */
  virtual void
  remove_column(const unsigned int k = 0);

  virtual void
  multiply_with_Q(VectorType &y, const Vector<Number> &x) const;

  virtual void
  multiply_with_QT(Vector<Number> &y, const VectorType &x) const;

  virtual void
  multiply_with_A(VectorType &y, const Vector<Number> &x) const;

  virtual void
  multiply_with_AT(Vector<Number> &y, const VectorType &x) const;

  /**
   * Connect a slot to implement a custom check of linear dependency
   * during addition of a column.
   *
   * Here, @p u is the last column of the to-be R matrix, @p rho
   * is its diagonal and @p col_norm_sqr is the square of the $l2$ norm of the column.
   * The function should return <code>true</code> if the new column is
   * linearly independent.
   */
  boost::signals2::connection
  connect_append_column_slot(
    const std::function<bool(const Vector<Number> &u,
                             const Number &        rho2,
                             const Number &        col_norm_sqr)> &slot);

private:
  /**
   * Apply givens rotation in the `(i,k)`-plane to zero out $R(k,k)$.
   */
  void
  apply_givens_rotation(const unsigned int i, const unsigned int k);

  /**
   * Signal used to decide if the new column is linear dependent.
   *
   * Here, @p u is the last column of the to-be R matrix, @p rho
   * is its diagonal and @p col_norm_sqr is the square of the $l2$ norm of the column.
   * The function should return <code>true</code> if the new column is
   * linearly independent.
   */
  boost::signals2::signal<bool(const Vector<Number> &u,
                               const Number &        rho,
                               const Number &        col_norm_sqr)>
    column_signal;
};

// -------------------  inline and template functions ----------------
#ifndef DOXYGEN

namespace internal
{
  namespace QRImplementation
  {
    // We want to avoid including our own LAPACK wrapper header in any external
    // headers to avoid possible conflicts with other packages that may define
    // their own such header. At the same time we want to be able to call some
    // LAPACK functions from the template functions below. To resolve both
    // problems define some extra wrappers here that can be in the header:
    template <typename Number>
    void
    call_trmv(const char            uplo,
              const char            trans,
              const char            diag,
              const types::blas_int n,
              const Number *        a,
              const types::blas_int lda,
              Number *              x,
              const types::blas_int incx);

    template <typename Number>
    void
    call_trtrs(const char            uplo,
               const char            trans,
               const char            diag,
               const types::blas_int n,
               const types::blas_int nrhs,
               const Number *        a,
               const types::blas_int lda,
               Number *              b,
               const types::blas_int ldb,
               types::blas_int *     info);
  } // namespace QRImplementation
} // namespace internal



template <typename VectorType>
BaseQR<VectorType>::BaseQR()
  : current_size(0)
{
  R.set_property(LAPACKSupport::upper_triangular);
}



template <typename VectorType>
unsigned int
BaseQR<VectorType>::size() const
{
  return current_size;
}



template <typename VectorType>
const LAPACKFullMatrix<typename BaseQR<VectorType>::Number> &
BaseQR<VectorType>::get_R() const
{
  return R;
}



template <typename VectorType>
void
BaseQR<VectorType>::solve(Vector<Number> &      x,
                          const Vector<Number> &y,
                          const bool            transpose) const
{
  Assert(x.size() == this->current_size,
         ExcDimensionMismatch(x.size(), this->current_size));
  Assert(y.size() == this->current_size,
         ExcDimensionMismatch(y.size(), this->current_size));

  // copy if the two vectors are not the same
  if (&x != &y)
    x = y;

  const int lda   = this->current_size;
  const int ldb   = this->current_size;
  const int N     = this->current_size;
  const int n_rhs = 1;
  int       info  = 0;
  internal::QRImplementation::call_trtrs('U',
                                         transpose ? 'T' : 'N',
                                         'N',
                                         N,
                                         n_rhs,
                                         &this->R(0, 0),
                                         lda,
                                         &x(0),
                                         ldb,
                                         &info);
}



template <typename VectorType>
void
BaseQR<VectorType>::multiply_with_cols(VectorType &          y,
                                       const Vector<Number> &x) const
{
  Assert(x.size() == this->current_size,
         ExcDimensionMismatch(x.size(), this->current_size));

  y = 0.;
  for (unsigned int j = 0; j < this->current_size; ++j)
    y.add(x[j], *this->columns[j]);
}



template <typename VectorType>
void
BaseQR<VectorType>::multiply_with_colsT(Vector<Number> &  y,
                                        const VectorType &x) const
{
  Assert(y.size() == this->current_size,
         ExcDimensionMismatch(y.size(), this->current_size));

  for (unsigned int j = 0; j < this->current_size; ++j)
    y[j] = (*this->columns[j]) * x;
}



template <class VectorType>
boost::signals2::connection
BaseQR<VectorType>::connect_givens_slot(
  const std::function<void(const unsigned int i,
                           const unsigned int j,
                           const std::array<Number, 3> &)> &slot)
{
  return givens_signal.connect(slot);
}



template <class VectorType>
boost::signals2::connection
ImplicitQR<VectorType>::connect_append_column_slot(
  const std::function<bool(const Vector<Number> &u,
                           const Number &        rho,
                           const Number &        col_norm_sqr)> &slot)
{
  return column_signal.connect(slot);
}



template <typename VectorType>
ImplicitQR<VectorType>::ImplicitQR()
  : BaseQR<VectorType>()
{}



template <typename VectorType>
bool
ImplicitQR<VectorType>::append_column(const VectorType &column)
{
  if (this->current_size == 0)
    {
      this->R.grow_or_shrink(this->current_size + 1);
      this->columns.push_back(std::make_unique<VectorType>(column));
      this->R(0, 0) = column.l2_norm();
      ++this->current_size;
    }
  else
    {
      // first get scalar products with A^T
      Vector<Number> u(this->current_size);
      this->multiply_with_AT(u, column);

      // now solve R^T x = (A^T * column)
      const int lda   = this->current_size;
      const int ldb   = this->current_size;
      const int N     = this->current_size;
      const int n_rhs = 1;
      int       info  = 0;
      internal::QRImplementation::call_trtrs(
        'U', 'T', 'N', N, n_rhs, &this->R(0, 0), lda, &u(0), ldb, &info);

      // finally get the diagonal element:
      // rho2 = |column|^2 - |u|^2
      const Number column_norm_sqr = column.norm_sqr();
      const Number rho2            = column_norm_sqr - u.norm_sqr();
      const bool   linearly_independent =
        column_signal.empty() ? rho2 > 0 :
                                column_signal(u, rho2, column_norm_sqr).get();

      // bail out if it turns out to be linearly dependent
      if (!linearly_independent)
        return false;

      // at this point we update is successful and we can enlarge R
      // and store the column:
      this->columns.push_back(std::make_unique<VectorType>(column));
      this->R.grow_or_shrink(this->current_size + 1);
      this->R(this->current_size, this->current_size) = std::sqrt(rho2);
      for (unsigned int i = 0; i < this->current_size; ++i)
        this->R(i, this->current_size) = u(i);

      this->current_size++;
    }

  return true;
}



template <typename VectorType>
void
ImplicitQR<VectorType>::apply_givens_rotation(const unsigned int i,
                                              const unsigned int k)
{
  AssertIndexRange(i, k);
  AssertIndexRange(k, this->current_size);
  const std::array<Number, 3> csr =
    dealii::Utilities::LinearAlgebra::givens_rotation<Number>(this->R(i, k),
                                                              this->R(k, k));

  // first, set k'th column:
  this->R(i, k) = csr[2];
  this->R(k, k) = 0.;
  // now do the rest:
  for (unsigned int j = 0; j < this->R.n(); ++j)
    if (j != k)
      {
        const Number t = this->R(i, j);
        this->R(i, j)  = csr[0] * this->R(i, j) + csr[1] * this->R(k, j);
        this->R(k, j)  = -csr[1] * t + csr[0] * this->R(k, j);
      }

  if (!this->givens_signal.empty())
    this->givens_signal(i, k, csr);
}



template <typename VectorType>
void
ImplicitQR<VectorType>::remove_column(const unsigned int k)
{
  // before actually removing a column from Q and resizing R,
  // apply givens rotations to bring H into upper triangular form:
  for (unsigned int j = k + 1; j < this->R.n(); ++j)
    {
      const unsigned int i = j - 1;
      apply_givens_rotation(i, j);
    }

  // remove last row and k-th column
  --this->current_size;
  this->R.remove_row_and_column(this->current_size, k);

  // Finally remove the column from A
  this->columns.erase(this->columns.begin() + k);
}



template <typename VectorType>
void
ImplicitQR<VectorType>::multiply_with_Q(VectorType &          y,
                                        const Vector<Number> &x) const
{
  // A = QR
  // A R^{-1} = Q
  Vector<Number> x1 = x;
  BaseQR<VectorType>::solve(x1, x1, false);
  multiply_with_A(y, x1);
}



template <typename VectorType>
void
ImplicitQR<VectorType>::multiply_with_QT(Vector<Number> &  y,
                                         const VectorType &x) const
{
  // A = QR
  // A^T = R^T Q^T
  // {R^T}^{-1} A^T = Q^T
  multiply_with_AT(y, x);
  BaseQR<VectorType>::solve(y, y, true);
}



template <typename VectorType>
void
ImplicitQR<VectorType>::multiply_with_A(VectorType &          y,
                                        const Vector<Number> &x) const
{
  BaseQR<VectorType>::multiply_with_cols(y, x);
}



template <typename VectorType>
void
ImplicitQR<VectorType>::multiply_with_AT(Vector<Number> &  y,
                                         const VectorType &x) const
{
  BaseQR<VectorType>::multiply_with_colsT(y, x);
}



template <typename VectorType>
QR<VectorType>::QR()
  : BaseQR<VectorType>()
{}



template <typename VectorType>
bool
QR<VectorType>::append_column(const VectorType &column)
{
  // resize R:
  this->R.grow_or_shrink(this->current_size + 1);
  this->columns.push_back(std::make_unique<VectorType>(column));

  // now a Gram-Schmidt part: orthonormalize the new column
  // against everything we have so far:
  auto &last_col = *this->columns.back();
  for (unsigned int i = 0; i < this->current_size; ++i)
    {
      const auto &i_col              = *this->columns[i];
      this->R(i, this->current_size) = i_col * last_col;
      last_col.add(-this->R(i, this->current_size), i_col);
    }

  this->R(this->current_size, this->current_size) = last_col.l2_norm();

  Assert(this->R(this->current_size, this->current_size) > 0.,
         ExcDivideByZero());
  last_col *= 1. / this->R(this->current_size, this->current_size);

  ++this->current_size;
  return true;
}



template <typename VectorType>
void
QR<VectorType>::apply_givens_rotation(const unsigned int i,
                                      const unsigned int k)
{
  AssertIndexRange(i, k);
  AssertIndexRange(k, this->current_size);
  const std::array<Number, 3> csr =
    dealii::Utilities::LinearAlgebra::givens_rotation<Number>(this->R(i, k),
                                                              this->R(k, k));

  // first, set k'th column:
  this->R(i, k) = csr[2];
  this->R(k, k) = 0.;
  // now do the rest:
  for (unsigned int j = 0; j < this->R.n(); ++j)
    if (j != k)
      {
        const Number t = this->R(i, j);
        this->R(i, j)  = csr[0] * this->R(i, j) + csr[1] * this->R(k, j);
        this->R(k, j)  = -csr[1] * t + csr[0] * this->R(k, j);
      }

  // now adjust i,k columns due to multiplication with the
  // transpose Givens matrix from right:
  auto &col_i = *this->columns[i];
  auto &col_k = *this->columns[k];
  // save column i:
  tmp = col_i;
  col_i.sadd(csr[0], csr[1], col_k);
  col_k.sadd(csr[0], -csr[1], tmp);

  if (!this->givens_signal.empty())
    this->givens_signal(i, k, csr);
}



template <typename VectorType>
void
QR<VectorType>::remove_column(const unsigned int k)
{
  AssertIndexRange(k, this->current_size);
  Assert(this->current_size > 0,
         ExcMessage("Can not remove a column if QR is empty"));
  // apply a sequence of Givens rotations
  // see section 6.5 "Updating matrix factorizations" in Golub 2013, Matrix
  // computations

  // So we want to have QR for \tilde A \in R^{m*(n-1)}
  // if we remove the column k, we end up with upper Hessenberg matrix
  //      x x x x x
  //        x x x x
  // H =      x x x
  //          x x x
  //            x x
  //              x
  // where k = 2 (3rd column), m = 7, n = 6
  //
  // before actually removing a column from Q and resizing R,
  // apply givens rotations to bring H into upper triangular form:
  for (unsigned int j = k + 1; j < this->R.n(); ++j)
    {
      const unsigned int i = j - 1;
      apply_givens_rotation(i, j);
    }

  // now we can throw away the column from Q and adjust R
  // since we do thin-QR, after Givens rotations we need to throw
  // away the last column:
  const unsigned int size_minus_1 = this->columns.size() - 1;
  this->columns.erase(this->columns.begin() + size_minus_1);

  // remove last row and k-th column
  --this->current_size;
  this->R.remove_row_and_column(this->current_size, k);
}



template <typename VectorType>
void
QR<VectorType>::multiply_with_Q(VectorType &y, const Vector<Number> &x) const
{
  BaseQR<VectorType>::multiply_with_cols(y, x);
}



template <typename VectorType>
void
QR<VectorType>::multiply_with_QT(Vector<Number> &y, const VectorType &x) const
{
  BaseQR<VectorType>::multiply_with_colsT(y, x);
}



template <typename VectorType>
void
QR<VectorType>::multiply_with_A(VectorType &y, const Vector<Number> &x) const
{
  Vector<Number> x1   = x;
  const int      N    = this->current_size;
  const int      lda  = N;
  const int      incx = 1;
  internal::QRImplementation::call_trmv(
    'U', 'N', 'N', N, &this->R(0, 0), lda, &x1[0], incx);

  multiply_with_Q(y, x1);
}



template <typename VectorType>
void
QR<VectorType>::multiply_with_AT(Vector<Number> &y, const VectorType &x) const
{
  multiply_with_QT(y, x);

  const int N    = this->current_size;
  const int lda  = N;
  const int incx = 1;
  internal::QRImplementation::call_trmv(
    'U', 'T', 'N', N, &this->R(0, 0), lda, &y[0], incx);
}

#endif // no DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
