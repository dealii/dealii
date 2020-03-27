// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2019 by the deal.II authors
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

#include <deal.II/base/numbers.h>

#include <deal.II/lac/blas_extension_templates.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/utilities.h>
#include <deal.II/lac/vector.h>

#include <iomanip>
#include <iostream>

DEAL_II_NAMESPACE_OPEN

using namespace LAPACKSupport;

namespace internal
{
  namespace LAPACKFullMatrixImplementation
  {
    // ZGEEV/CGEEV and DGEEV/SGEEV need different work arrays and different
    // output arrays for eigenvalues. This makes working with generic scalar
    // types a bit difficult. To get around this, geev_helper has the same
    // signature for real and complex arguments, but it ignores some
    // parameters when called with a real type and ignores different
    // parameters when called with a complex type.
    template <typename T>
    void
    geev_helper(const char            vl,
                const char            vr,
                AlignedVector<T> &    matrix,
                const types::blas_int n_rows,
                std::vector<T> &      real_part_eigenvalues,
                std::vector<T> &      imag_part_eigenvalues,
                std::vector<T> &      left_eigenvectors,
                std::vector<T> &      right_eigenvectors,
                std::vector<T> &      real_work,
                std::vector<T> & /*complex_work*/,
                const types::blas_int work_flag,
                types::blas_int &     info)
    {
      static_assert(std::is_same<T, double>::value ||
                      std::is_same<T, float>::value,
                    "Only implemented for double and float");
      Assert(matrix.size() == static_cast<std::size_t>(n_rows * n_rows),
             ExcInternalError());
      Assert(static_cast<std::size_t>(n_rows) <= real_part_eigenvalues.size(),
             ExcInternalError());
      Assert(static_cast<std::size_t>(n_rows) <= imag_part_eigenvalues.size(),
             ExcInternalError());
      if (vl == 'V')
        Assert(static_cast<std::size_t>(n_rows * n_rows) <=
                 left_eigenvectors.size(),
               ExcInternalError());
      if (vr == 'V')
        Assert(static_cast<std::size_t>(n_rows * n_rows) <=
                 right_eigenvectors.size(),
               ExcInternalError());
      Assert(work_flag == -1 ||
               static_cast<std::size_t>(2 * n_rows) <= real_work.size(),
             ExcInternalError());
      Assert(work_flag == -1 || std::max<long int>(1, 3 * n_rows) <= work_flag,
             ExcInternalError());
      geev(&vl,
           &vr,
           &n_rows,
           matrix.data(),
           &n_rows,
           real_part_eigenvalues.data(),
           imag_part_eigenvalues.data(),
           left_eigenvectors.data(),
           &n_rows,
           right_eigenvectors.data(),
           &n_rows,
           real_work.data(),
           &work_flag,
           &info);
    }

    template <typename T>
    void
    geev_helper(const char                      vl,
                const char                      vr,
                AlignedVector<std::complex<T>> &matrix,
                const types::blas_int           n_rows,
                std::vector<T> & /*real_part_eigenvalues*/,
                std::vector<std::complex<T>> &eigenvalues,
                std::vector<std::complex<T>> &left_eigenvectors,
                std::vector<std::complex<T>> &right_eigenvectors,
                std::vector<std::complex<T>> &complex_work,
                std::vector<T> &              real_work,
                const types::blas_int         work_flag,
                types::blas_int &             info)
    {
      static_assert(
        std::is_same<T, double>::value || std::is_same<T, float>::value,
        "Only implemented for std::complex<double> and std::complex<float>");
      Assert(matrix.size() == static_cast<std::size_t>(n_rows * n_rows),
             ExcInternalError());
      Assert(static_cast<std::size_t>(n_rows) <= eigenvalues.size(),
             ExcInternalError());
      if (vl == 'V')
        Assert(static_cast<std::size_t>(n_rows * n_rows) <=
                 left_eigenvectors.size(),
               ExcInternalError());
      if (vr == 'V')
        Assert(static_cast<std::size_t>(n_rows * n_rows) <=
                 right_eigenvectors.size(),
               ExcInternalError());
      Assert(std::max<std::size_t>(1, work_flag) <= real_work.size(),
             ExcInternalError());
      Assert(work_flag == -1 ||
               std::max<long int>(1, 2 * n_rows) <= (work_flag),
             ExcInternalError());

      geev(&vl,
           &vr,
           &n_rows,
           matrix.data(),
           &n_rows,
           eigenvalues.data(),
           left_eigenvectors.data(),
           &n_rows,
           right_eigenvectors.data(),
           &n_rows,
           complex_work.data(),
           &work_flag,
           real_work.data(),
           &info);
    }



    template <typename T>
    void
    gesdd_helper(const char            job,
                 const types::blas_int n_rows,
                 const types::blas_int n_cols,
                 AlignedVector<T> &    matrix,
                 std::vector<T> &      singular_values,
                 AlignedVector<T> &    left_vectors,
                 AlignedVector<T> &    right_vectors,
                 std::vector<T> &      real_work,
                 std::vector<T> & /*complex work*/,
                 std::vector<types::blas_int> &integer_work,
                 const types::blas_int         work_flag,
                 types::blas_int &             info)
    {
      Assert(job == 'A' || job == 'S' || job == 'O' || job == 'N',
             ExcInternalError());
      Assert(static_cast<std::size_t>(n_rows * n_cols) == matrix.size(),
             ExcInternalError());
      Assert(std::min<std::size_t>(n_rows, n_cols) <= singular_values.size(),
             ExcInternalError());
      Assert(8 * std::min<std::size_t>(n_rows, n_cols) <= integer_work.size(),
             ExcInternalError());
      Assert(work_flag == -1 ||
               static_cast<std::size_t>(work_flag) <= real_work.size(),
             ExcInternalError());
      gesdd(&job,
            &n_rows,
            &n_cols,
            matrix.data(),
            &n_rows,
            singular_values.data(),
            left_vectors.data(),
            &n_rows,
            right_vectors.data(),
            &n_cols,
            real_work.data(),
            &work_flag,
            integer_work.data(),
            &info);
    }



    template <typename T>
    void
    gesdd_helper(const char                      job,
                 const types::blas_int           n_rows,
                 const types::blas_int           n_cols,
                 AlignedVector<std::complex<T>> &matrix,
                 std::vector<T> &                singular_values,
                 AlignedVector<std::complex<T>> &left_vectors,
                 AlignedVector<std::complex<T>> &right_vectors,
                 std::vector<std::complex<T>> &  work,
                 std::vector<T> &                real_work,
                 std::vector<types::blas_int> &  integer_work,
                 const types::blas_int &         work_flag,
                 types::blas_int &               info)
    {
      Assert(job == 'A' || job == 'S' || job == 'O' || job == 'N',
             ExcInternalError());
      Assert(static_cast<std::size_t>(n_rows * n_cols) == matrix.size(),
             ExcInternalError());
      Assert(static_cast<std::size_t>(std::min(n_rows, n_cols)) <=
               singular_values.size(),
             ExcInternalError());
      Assert(8 * std::min<std::size_t>(n_rows, n_cols) <= integer_work.size(),
             ExcInternalError());
      Assert(work_flag == -1 ||
               static_cast<std::size_t>(work_flag) <= real_work.size(),
             ExcInternalError());

      gesdd(&job,
            &n_rows,
            &n_cols,
            matrix.data(),
            &n_rows,
            singular_values.data(),
            left_vectors.data(),
            &n_rows,
            right_vectors.data(),
            &n_cols,
            work.data(),
            &work_flag,
            real_work.data(),
            integer_work.data(),
            &info);
    }
  } // namespace LAPACKFullMatrixImplementation
} // namespace internal

template <typename number>
LAPACKFullMatrix<number>::LAPACKFullMatrix(const size_type n)
  : TransposeTable<number>(n, n)
  , state(matrix)
  , property(general)
{}


template <typename number>
LAPACKFullMatrix<number>::LAPACKFullMatrix(const size_type m, const size_type n)
  : TransposeTable<number>(m, n)
  , state(matrix)
  , property(general)
{}


template <typename number>
LAPACKFullMatrix<number>::LAPACKFullMatrix(const LAPACKFullMatrix &M)
  : TransposeTable<number>(M)
  , state(matrix)
  , property(general)
{}


template <typename number>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator=(const LAPACKFullMatrix<number> &M)
{
  TransposeTable<number>::operator=(M);
  state                           = M.state;
  property                        = M.property;
  return *this;
}



template <typename number>
void
LAPACKFullMatrix<number>::reinit(const size_type n)
{
  this->TransposeTable<number>::reinit(n, n);
  state = LAPACKSupport::matrix;
}



template <typename number>
void
LAPACKFullMatrix<number>::grow_or_shrink(const size_type n)
{
  const size_type        s = std::min(std::min(this->m(), n), this->n());
  TransposeTable<number> copy(std::move(*this));
  this->TransposeTable<number>::reinit(n, n);
  for (size_type i = 0; i < s; ++i)
    for (size_type j = 0; j < s; ++j)
      (*this)(i, j) = copy(i, j);
}



template <typename number>
void
LAPACKFullMatrix<number>::apply_givens_rotation(
  const std::array<number, 3> &csr,
  const size_type              i,
  const size_type              k,
  const bool                   left)
{
  auto &A = *this;
  // see Golub 2013 "Matrix computations", p241 5.1.9 Applying Givens
  // Rotations but note the difference in notation, namely the sign of s: we
  // have G * A, where G[1,1] = s
  if (left)
    {
      for (size_type j = 0; j < A.n(); ++j)
        {
          const number t = A(i, j);
          A(i, j)        = csr[0] * A(i, j) + csr[1] * A(k, j);
          A(k, j)        = -csr[1] * t + csr[0] * A(k, j);
        }
    }
  else
    {
      for (size_type j = 0; j < A.m(); ++j)
        {
          const number t = A(j, i);
          A(j, i)        = csr[0] * A(j, i) + csr[1] * A(j, k);
          A(j, k)        = -csr[1] * t + csr[0] * A(j, k);
        }
    }
}



template <typename number>
void
LAPACKFullMatrix<number>::remove_row_and_column(const size_type row,
                                                const size_type col)
{
  AssertIndexRange(row, this->m());
  AssertIndexRange(col, this->n());

  const size_type nrows = this->m() - 1;
  const size_type ncols = this->n() - 1;

  TransposeTable<number> copy(std::move(*this));
  this->TransposeTable<number>::reinit(nrows, ncols);

  for (size_type j = 0; j < ncols; ++j)
    {
      const size_type jj = (j < col ? j : j + 1);
      for (size_type i = 0; i < nrows; ++i)
        {
          const size_type ii = (i < row ? i : i + 1);
          (*this)(i, j)      = copy(ii, jj);
        }
    }
}



template <typename number>
void
LAPACKFullMatrix<number>::reinit(const size_type m, const size_type n)
{
  this->TransposeTable<number>::reinit(m, n);
  state = LAPACKSupport::matrix;
}


template <typename number>
template <typename number2>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator=(const FullMatrix<number2> &M)
{
  Assert(this->m() == M.m(), ExcDimensionMismatch(this->m(), M.m()));
  Assert(this->n() == M.n(), ExcDimensionMismatch(this->n(), M.n()));
  for (size_type i = 0; i < this->m(); ++i)
    for (size_type j = 0; j < this->n(); ++j)
      (*this)(i, j) = M(i, j);

  state    = LAPACKSupport::matrix;
  property = LAPACKSupport::general;
  return *this;
}


template <typename number>
template <typename number2>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator=(const SparseMatrix<number2> &M)
{
  Assert(this->m() == M.n(), ExcDimensionMismatch(this->m(), M.n()));
  Assert(this->n() == M.m(), ExcDimensionMismatch(this->n(), M.m()));
  for (size_type i = 0; i < this->m(); ++i)
    for (size_type j = 0; j < this->n(); ++j)
      (*this)(i, j) = M.el(i, j);

  state    = LAPACKSupport::matrix;
  property = LAPACKSupport::general;
  return *this;
}


template <typename number>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator=(const number d)
{
  (void)d;
  Assert(d == number(0), ExcScalarAssignmentOnlyForZeroValue());

  if (this->n_elements() != 0)
    this->reset_values();

  state = LAPACKSupport::matrix;
  return *this;
}


template <typename number>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator*=(const number factor)
{
  Assert(state == LAPACKSupport::matrix ||
           state == LAPACKSupport::inverse_matrix,
         ExcState(state));

  AssertIsFinite(factor);
  const char            type  = 'G';
  const number          cfrom = 1.;
  const types::blas_int m     = this->m();
  const types::blas_int n     = this->n();
  const types::blas_int lda   = this->m();
  types::blas_int       info  = 0;
  // kl and ku will not be referenced for type = G (dense matrices).
  const types::blas_int kl     = 0;
  number *              values = this->values.data();

  lascl(&type, &kl, &kl, &cfrom, &factor, &m, &n, values, &lda, &info);

  // Negative return value implies a wrong argument. This should be internal.
  Assert(info >= 0, ExcInternalError());

  return *this;
}


template <typename number>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator/=(const number factor)
{
  Assert(state == LAPACKSupport::matrix ||
           state == LAPACKSupport::inverse_matrix,
         ExcState(state));

  AssertIsFinite(factor);
  Assert(factor != number(0.), ExcZero());

  const char            type = 'G';
  const number          cto  = 1.;
  const types::blas_int m    = this->m();
  const types::blas_int n    = this->n();
  const types::blas_int lda  = this->m();
  types::blas_int       info = 0;
  // kl and ku will not be referenced for type = G (dense matrices).
  const types::blas_int kl     = 0;
  number *              values = this->values.data();

  lascl(&type, &kl, &kl, &factor, &cto, &m, &n, values, &lda, &info);

  // Negative return value implies a wrong argument. This should be internal.
  Assert(info >= 0, ExcInternalError());

  return *this;
}



template <typename number>
void
LAPACKFullMatrix<number>::add(const number a, const LAPACKFullMatrix<number> &A)
{
  Assert(state == LAPACKSupport::matrix ||
           state == LAPACKSupport::inverse_matrix,
         ExcState(state));

  Assert(m() == A.m(), ExcDimensionMismatch(m(), A.m()));
  Assert(n() == A.n(), ExcDimensionMismatch(n(), A.n()));

  AssertIsFinite(a);

  // BLAS does not offer functions to add matrices.
  // LapackFullMatrix is stored in contiguous array
  // ==> use BLAS 1 for adding vectors
  const types::blas_int n        = this->m() * this->n();
  const types::blas_int inc      = 1;
  number *              values   = this->values.data();
  const number *        values_A = A.values.data();

  axpy(&n, &a, values_A, &inc, values, &inc);
}



namespace
{
  template <typename number>
  void
  cholesky_rank1(LAPACKFullMatrix<number> &A,
                 const number              a,
                 const Vector<number> &    v)
  {
    const typename LAPACKFullMatrix<number>::size_type N = A.n();
    Vector<number>                                     z(v);
    // Cholesky update / downdate, see
    // 6.5.4 Cholesky Updating and Downdating, Golub 2013 Matrix computations
    // Note that potrf() is called with LAPACKSupport::L , so the
    // factorization is stored in lower triangular part.
    // Also see discussion here
    // http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=2646
    if (a > 0.)
      {
        // simple update via a sequence of Givens rotations.
        // Observe that
        //
        //       | L^T |T  | L^T |
        // A <-- |     |   |     | = L L^T + z z^T
        //       | z^T |   | z^T |
        //
        // so we can get updated factor by doing a sequence of Givens
        // rotations to make the matrix lower-triangular
        // Also see LINPACK's dchud http://www.netlib.org/linpack/dchud.f
        z *= std::sqrt(a);
        for (typename LAPACKFullMatrix<number>::size_type k = 0; k < N; ++k)
          {
            const std::array<number, 3> csr =
              Utilities::LinearAlgebra::givens_rotation(A(k, k), z(k));
            A(k, k) = csr[2];
            for (typename LAPACKFullMatrix<number>::size_type i = k + 1; i < N;
                 ++i)
              {
                const number t = A(i, k);
                A(i, k)        = csr[0] * A(i, k) + csr[1] * z(i);
                z(i)           = -csr[1] * t + csr[0] * z(i);
              }
          }
      }
    else
      {
        // downdating is not always possible as we may end up with
        // negative definite matrix. If it's possible, then it boils
        // down to application of hyperbolic rotations.
        // Observe that
        //
        //       | L^T |T      | L^T |
        // A <-- |     |   S   |     | = L L^T - z z^T
        //       | z^T |       | z^T |
        //
        //       |In  0 |
        // S :=  |      |
        //       |0  -1 |
        //
        // We are looking for H which is S-orthogonal (HSH^T=S) and
        // can restore upper-triangular factor of the factorization of A above.
        // We will use Hyperbolic rotations to do the job
        //
        // | c  -s |   | x1 |   | r |
        // |       | = |    | = |   |,   c^2 - s^2 = 1
        // |-s   c |   | x2 |   | 0 |
        //
        // which have real solution only if x2 <= x1.
        // See also Linpack's http://www.netlib.org/linpack/dchdd.f and
        // https://infoscience.epfl.ch/record/161468/files/cholupdate.pdf and
        // "Analysis of a recursive Least Squares Hyperbolic Rotation Algorithm
        // for Signal Processing", Alexander, Pan, Plemmons, 1988.
        z *= std::sqrt(-a);
        for (typename LAPACKFullMatrix<number>::size_type k = 0; k < N; ++k)
          {
            const std::array<number, 3> csr =
              Utilities::LinearAlgebra::hyperbolic_rotation(A(k, k), z(k));
            A(k, k) = csr[2];
            for (typename LAPACKFullMatrix<number>::size_type i = k + 1; i < N;
                 ++i)
              {
                const number t = A(i, k);
                A(i, k)        = csr[0] * A(i, k) - csr[1] * z(i);
                z(i)           = -csr[1] * t + csr[0] * z(i);
              }
          }
      }
  }


  template <typename number>
  void
  cholesky_rank1(LAPACKFullMatrix<std::complex<number>> & /*A*/,
                 const std::complex<number> /*a*/,
                 const Vector<std::complex<number>> & /*v*/)
  {
    AssertThrow(false, ExcNotImplemented());
  }
} // namespace



template <typename number>
void
LAPACKFullMatrix<number>::rank1_update(const number a, const Vector<number> &v)
{
  Assert(property == LAPACKSupport::symmetric, ExcProperty(property));

  Assert(n() == m(), ExcInternalError());
  Assert(m() == v.size(), ExcDimensionMismatch(m(), v.size()));

  AssertIsFinite(a);

  if (state == LAPACKSupport::matrix)
    {
      {
        const types::blas_int N    = this->m();
        const char            uplo = LAPACKSupport::U;
        const types::blas_int lda  = N;
        const types::blas_int incx = 1;

        syr(&uplo, &N, &a, v.begin(), &incx, this->values.begin(), &lda);
      }

      const size_type N = this->m();
      // FIXME: we should really only update upper or lower triangular parts
      // of symmetric matrices and make sure the interface is consistent,
      // for example operator(i,j) gives correct results regardless of storage.
      for (size_type i = 0; i < N; ++i)
        for (size_type j = 0; j < i; ++j)
          (*this)(i, j) = (*this)(j, i);
    }
  else if (state == LAPACKSupport::cholesky)
    {
      cholesky_rank1(*this, a, v);
    }
  else
    AssertThrow(false, ExcState(state));
}



template <typename number>
void
LAPACKFullMatrix<number>::vmult(Vector<number> &      w,
                                const Vector<number> &v,
                                const bool            adding) const
{
  const types::blas_int mm    = this->m();
  const types::blas_int nn    = this->n();
  const number          alpha = 1.;
  const number          beta  = (adding ? 1. : 0.);
  const number          null  = 0.;

  // use trmv for triangular matrices
  if ((property == upper_triangular || property == lower_triangular) &&
      (mm == nn) && state == matrix)
    {
      Assert(adding == false, ExcNotImplemented());

      AssertDimension(v.size(), this->n());
      AssertDimension(w.size(), this->m());

      const char diag  = 'N';
      const char trans = 'N';
      const char uplo =
        (property == upper_triangular ? LAPACKSupport::U : LAPACKSupport::L);

      w = v;

      const types::blas_int N    = mm;
      const types::blas_int lda  = N;
      const types::blas_int incx = 1;

      trmv(
        &uplo, &trans, &diag, &N, this->values.data(), &lda, w.data(), &incx);

      return;
    }

  switch (state)
    {
      case matrix:
      case inverse_matrix:
        {
          AssertDimension(v.size(), this->n());
          AssertDimension(w.size(), this->m());

          gemv("N",
               &mm,
               &nn,
               &alpha,
               this->values.data(),
               &mm,
               v.data(),
               &one,
               &beta,
               w.data(),
               &one);
          break;
        }
      case svd:
        {
          std::lock_guard<std::mutex> lock(mutex);
          AssertDimension(v.size(), this->n());
          AssertDimension(w.size(), this->m());
          // Compute V^T v
          work.resize(std::max(mm, nn));
          gemv("N",
               &nn,
               &nn,
               &alpha,
               svd_vt->values.data(),
               &nn,
               v.data(),
               &one,
               &null,
               work.data(),
               &one);
          // Multiply by singular values
          for (size_type i = 0; i < wr.size(); ++i)
            work[i] *= wr[i];
          // Multiply with U
          gemv("N",
               &mm,
               &mm,
               &alpha,
               svd_u->values.data(),
               &mm,
               work.data(),
               &one,
               &beta,
               w.data(),
               &one);
          break;
        }
      case inverse_svd:
        {
          std::lock_guard<std::mutex> lock(mutex);
          AssertDimension(w.size(), this->n());
          AssertDimension(v.size(), this->m());
          // Compute U^T v
          work.resize(std::max(mm, nn));
          gemv("T",
               &mm,
               &mm,
               &alpha,
               svd_u->values.data(),
               &mm,
               v.data(),
               &one,
               &null,
               work.data(),
               &one);
          // Multiply by singular values
          for (size_type i = 0; i < wr.size(); ++i)
            work[i] *= wr[i];
          // Multiply with V
          gemv("T",
               &nn,
               &nn,
               &alpha,
               svd_vt->values.data(),
               &nn,
               work.data(),
               &one,
               &beta,
               w.data(),
               &one);
          break;
        }
      default:
        Assert(false, ExcState(state));
    }
}


template <typename number>
void
LAPACKFullMatrix<number>::Tvmult(Vector<number> &      w,
                                 const Vector<number> &v,
                                 const bool            adding) const
{
  const types::blas_int mm    = this->m();
  const types::blas_int nn    = this->n();
  const number          alpha = 1.;
  const number          beta  = (adding ? 1. : 0.);
  const number          null  = 0.;

  // use trmv for triangular matrices
  if ((property == upper_triangular || property == lower_triangular) &&
      (mm == nn) && state == matrix)
    {
      Assert(adding == false, ExcNotImplemented());

      AssertDimension(v.size(), this->n());
      AssertDimension(w.size(), this->m());

      const char diag  = 'N';
      const char trans = 'T';
      const char uplo =
        (property == upper_triangular ? LAPACKSupport::U : LAPACKSupport::L);

      w = v;

      const types::blas_int N    = mm;
      const types::blas_int lda  = N;
      const types::blas_int incx = 1;

      trmv(
        &uplo, &trans, &diag, &N, this->values.data(), &lda, w.data(), &incx);

      return;
    }


  switch (state)
    {
      case matrix:
      case inverse_matrix:
        {
          AssertDimension(w.size(), this->n());
          AssertDimension(v.size(), this->m());

          gemv("T",
               &mm,
               &nn,
               &alpha,
               this->values.data(),
               &mm,
               v.data(),
               &one,
               &beta,
               w.data(),
               &one);
          break;
        }
      case svd:
        {
          std::lock_guard<std::mutex> lock(mutex);
          AssertDimension(w.size(), this->n());
          AssertDimension(v.size(), this->m());

          // Compute U^T v
          work.resize(std::max(mm, nn));
          gemv("T",
               &mm,
               &mm,
               &alpha,
               svd_u->values.data(),
               &mm,
               v.data(),
               &one,
               &null,
               work.data(),
               &one);
          // Multiply by singular values
          for (size_type i = 0; i < wr.size(); ++i)
            work[i] *= wr[i];
          // Multiply with V
          gemv("T",
               &nn,
               &nn,
               &alpha,
               svd_vt->values.data(),
               &nn,
               work.data(),
               &one,
               &beta,
               w.data(),
               &one);
          break;
        }
      case inverse_svd:
        {
          std::lock_guard<std::mutex> lock(mutex);
          AssertDimension(v.size(), this->n());
          AssertDimension(w.size(), this->m());

          // Compute V^T v
          work.resize(std::max(mm, nn));
          gemv("N",
               &nn,
               &nn,
               &alpha,
               svd_vt->values.data(),
               &nn,
               v.data(),
               &one,
               &null,
               work.data(),
               &one);
          // Multiply by singular values
          for (size_type i = 0; i < wr.size(); ++i)
            work[i] *= wr[i];
          // Multiply with U
          gemv("N",
               &mm,
               &mm,
               &alpha,
               svd_u->values.data(),
               &mm,
               work.data(),
               &one,
               &beta,
               w.data(),
               &one);
          break;
        }
      default:
        Assert(false, ExcState(state));
    }
}


template <typename number>
void
LAPACKFullMatrix<number>::vmult_add(Vector<number> &      w,
                                    const Vector<number> &v) const
{
  vmult(w, v, true);
}


template <typename number>
void
LAPACKFullMatrix<number>::Tvmult_add(Vector<number> &      w,
                                     const Vector<number> &v) const
{
  Tvmult(w, v, true);
}


template <typename number>
void
LAPACKFullMatrix<number>::mmult(LAPACKFullMatrix<number> &      C,
                                const LAPACKFullMatrix<number> &B,
                                const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(B.state));
  Assert(C.state == matrix || C.state == inverse_matrix, ExcState(C.state));
  Assert(this->n() == B.m(), ExcDimensionMismatch(this->n(), B.m()));
  Assert(C.n() == B.n(), ExcDimensionMismatch(C.n(), B.n()));
  Assert(C.m() == this->m(), ExcDimensionMismatch(this->m(), C.m()));
  const types::blas_int mm    = this->m();
  const types::blas_int nn    = B.n();
  const types::blas_int kk    = this->n();
  const number          alpha = 1.;
  const number          beta  = (adding ? 1. : 0.);

  gemm("N",
       "N",
       &mm,
       &nn,
       &kk,
       &alpha,
       this->values.data(),
       &mm,
       B.values.data(),
       &kk,
       &beta,
       C.values.data(),
       &mm);
}


template <typename number>
void
LAPACKFullMatrix<number>::mmult(FullMatrix<number> &            C,
                                const LAPACKFullMatrix<number> &B,
                                const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(B.state));
  Assert(this->n() == B.m(), ExcDimensionMismatch(this->n(), B.m()));
  Assert(C.n() == B.n(), ExcDimensionMismatch(C.n(), B.n()));
  Assert(C.m() == this->m(), ExcDimensionMismatch(this->m(), C.m()));
  const types::blas_int mm    = this->m();
  const types::blas_int nn    = B.n();
  const types::blas_int kk    = this->n();
  const number          alpha = 1.;
  const number          beta  = (adding ? 1. : 0.);

  // since FullMatrix stores the matrix in transposed order compared to this
  // matrix, compute B^T * A^T = (A * B)^T
  gemm("T",
       "T",
       &nn,
       &mm,
       &kk,
       &alpha,
       B.values.data(),
       &kk,
       this->values.data(),
       &mm,
       &beta,
       &C(0, 0),
       &nn);
}



template <typename number>
void
LAPACKFullMatrix<number>::Tmmult(LAPACKFullMatrix<number> &      C,
                                 const LAPACKFullMatrix<number> &B,
                                 const Vector<number> &          V,
                                 const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(B.state));
  Assert(C.state == matrix || C.state == inverse_matrix, ExcState(C.state));

  const LAPACKFullMatrix<number> &A = *this;

  Assert(A.m() == B.m(), ExcDimensionMismatch(A.m(), B.m()));
  Assert(C.n() == B.n(), ExcDimensionMismatch(C.n(), B.n()));
  Assert(C.m() == A.n(), ExcDimensionMismatch(A.n(), C.m()));
  Assert(V.size() == A.m(), ExcDimensionMismatch(A.m(), V.size()));

  const types::blas_int mm = A.n();
  const types::blas_int nn = B.n();
  const types::blas_int kk = B.m();

  // lapack does not have any triple product routines, including the case of
  // diagonal matrix in the middle, see
  // https://stackoverflow.com/questions/3548069/multiplying-three-matrices-in-blas-with-the-middle-one-being-diagonal
  // http://mathforum.org/kb/message.jspa?messageID=3546564

  std::lock_guard<std::mutex> lock(mutex);
  // First, get V*B into "work" array
  work.resize(kk * nn);
  // following http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=768#p2577
  // do left-multiplication manually. Note that Xscal would require to first
  // copy the input vector as multiplication is done inplace.
  for (types::blas_int j = 0; j < nn; ++j)
    for (types::blas_int i = 0; i < kk; ++i)
      {
        Assert(j * kk + i < static_cast<types::blas_int>(work.size()),
               ExcInternalError());
        work[j * kk + i] = V(i) * B(i, j);
      }

  // Now do the standard Tmmult:
  const number alpha = 1.;
  const number beta  = (adding ? 1. : 0.);

  gemm("T",
       "N",
       &mm,
       &nn,
       &kk,
       &alpha,
       this->values.data(),
       &kk,
       work.data(),
       &kk,
       &beta,
       C.values.data(),
       &mm);
}



template <typename number>
void
LAPACKFullMatrix<number>::transpose(LAPACKFullMatrix<number> &B) const
{
  const LAPACKFullMatrix<number> &A = *this;
  AssertDimension(A.m(), B.n());
  AssertDimension(A.n(), B.m());
  const types::blas_int m = B.m();
  const types::blas_int n = B.n();
#ifdef DEAL_II_LAPACK_WITH_MKL
  const number one = 1.;
  omatcopy('C', 'C', n, m, one, A.values.data(), n, B.values.data(), m);
#else
  for (types::blas_int i = 0; i < m; ++i)
    for (types::blas_int j = 0; j < n; ++j)
      B(i, j) = numbers::NumberTraits<number>::conjugate(A(j, i));
#endif
}


template <typename number>
void
LAPACKFullMatrix<number>::scale_rows(const Vector<number> &V)
{
  LAPACKFullMatrix<number> &A = *this;
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(V.size() == A.m(), ExcDimensionMismatch(A.m(), V.size()));

  const types::blas_int nn = A.n();
  const types::blas_int kk = A.m();
  for (types::blas_int j = 0; j < nn; ++j)
    for (types::blas_int i = 0; i < kk; ++i)
      A(i, j) *= V(i);
}



template <typename number>
void
LAPACKFullMatrix<number>::Tmmult(LAPACKFullMatrix<number> &      C,
                                 const LAPACKFullMatrix<number> &B,
                                 const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(B.state));
  Assert(C.state == matrix || C.state == inverse_matrix, ExcState(C.state));
  Assert(this->m() == B.m(), ExcDimensionMismatch(this->m(), B.m()));
  Assert(C.n() == B.n(), ExcDimensionMismatch(C.n(), B.n()));
  Assert(C.m() == this->n(), ExcDimensionMismatch(this->n(), C.m()));
  const types::blas_int mm    = this->n();
  const types::blas_int nn    = B.n();
  const types::blas_int kk    = B.m();
  const number          alpha = 1.;
  const number          beta  = (adding ? 1. : 0.);

  if (PointerComparison::equal(this, &B))
    {
      syrk(&LAPACKSupport::U,
           "T",
           &nn,
           &kk,
           &alpha,
           this->values.data(),
           &kk,
           &beta,
           C.values.data(),
           &nn);

      // fill-in lower triangular part
      for (types::blas_int j = 0; j < nn; ++j)
        for (types::blas_int i = 0; i < j; ++i)
          C(j, i) = C(i, j);

      C.property = symmetric;
    }
  else
    {
      gemm("T",
           "N",
           &mm,
           &nn,
           &kk,
           &alpha,
           this->values.data(),
           &kk,
           B.values.data(),
           &kk,
           &beta,
           C.values.data(),
           &mm);
    }
}


template <typename number>
void
LAPACKFullMatrix<number>::Tmmult(FullMatrix<number> &            C,
                                 const LAPACKFullMatrix<number> &B,
                                 const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(B.state));
  Assert(this->m() == B.m(), ExcDimensionMismatch(this->m(), B.m()));
  Assert(C.n() == B.n(), ExcDimensionMismatch(C.n(), B.n()));
  Assert(C.m() == this->n(), ExcDimensionMismatch(this->n(), C.m()));
  const types::blas_int mm    = this->n();
  const types::blas_int nn    = B.n();
  const types::blas_int kk    = B.m();
  const number          alpha = 1.;
  const number          beta  = (adding ? 1. : 0.);

  // since FullMatrix stores the matrix in transposed order compared to this
  // matrix, compute B^T * A = (A^T * B)^T
  gemm("T",
       "N",
       &nn,
       &mm,
       &kk,
       &alpha,
       B.values.data(),
       &kk,
       this->values.data(),
       &kk,
       &beta,
       &C(0, 0),
       &nn);
}


template <typename number>
void
LAPACKFullMatrix<number>::mTmult(LAPACKFullMatrix<number> &      C,
                                 const LAPACKFullMatrix<number> &B,
                                 const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(B.state));
  Assert(C.state == matrix || C.state == inverse_matrix, ExcState(C.state));
  Assert(this->n() == B.n(), ExcDimensionMismatch(this->n(), B.n()));
  Assert(C.n() == B.m(), ExcDimensionMismatch(C.n(), B.m()));
  Assert(C.m() == this->m(), ExcDimensionMismatch(this->m(), C.m()));
  const types::blas_int mm    = this->m();
  const types::blas_int nn    = B.m();
  const types::blas_int kk    = B.n();
  const number          alpha = 1.;
  const number          beta  = (adding ? 1. : 0.);

  if (PointerComparison::equal(this, &B))
    {
      syrk(&LAPACKSupport::U,
           "N",
           &nn,
           &kk,
           &alpha,
           this->values.data(),
           &nn,
           &beta,
           C.values.data(),
           &nn);

      // fill-in lower triangular part
      for (types::blas_int j = 0; j < nn; ++j)
        for (types::blas_int i = 0; i < j; ++i)
          C(j, i) = C(i, j);

      C.property = symmetric;
    }
  else
    {
      gemm("N",
           "T",
           &mm,
           &nn,
           &kk,
           &alpha,
           this->values.data(),
           &mm,
           B.values.data(),
           &nn,
           &beta,
           C.values.data(),
           &mm);
    }
}



template <typename number>
void
LAPACKFullMatrix<number>::mTmult(FullMatrix<number> &            C,
                                 const LAPACKFullMatrix<number> &B,
                                 const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(B.state));
  Assert(this->n() == B.n(), ExcDimensionMismatch(this->n(), B.n()));
  Assert(C.n() == B.m(), ExcDimensionMismatch(C.n(), B.m()));
  Assert(C.m() == this->m(), ExcDimensionMismatch(this->m(), C.m()));
  const types::blas_int mm    = this->m();
  const types::blas_int nn    = B.m();
  const types::blas_int kk    = B.n();
  const number          alpha = 1.;
  const number          beta  = (adding ? 1. : 0.);

  // since FullMatrix stores the matrix in transposed order compared to this
  // matrix, compute B * A^T = (A * B^T)^T
  gemm("N",
       "T",
       &nn,
       &mm,
       &kk,
       &alpha,
       B.values.data(),
       &nn,
       this->values.data(),
       &mm,
       &beta,
       &C(0, 0),
       &nn);
}


template <typename number>
void
LAPACKFullMatrix<number>::TmTmult(LAPACKFullMatrix<number> &      C,
                                  const LAPACKFullMatrix<number> &B,
                                  const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(B.state));
  Assert(C.state == matrix || C.state == inverse_matrix, ExcState(C.state));
  Assert(this->m() == B.n(), ExcDimensionMismatch(this->m(), B.n()));
  Assert(C.n() == B.m(), ExcDimensionMismatch(C.n(), B.m()));
  Assert(C.m() == this->n(), ExcDimensionMismatch(this->n(), C.m()));
  const types::blas_int mm    = this->n();
  const types::blas_int nn    = B.m();
  const types::blas_int kk    = B.n();
  const number          alpha = 1.;
  const number          beta  = (adding ? 1. : 0.);

  gemm("T",
       "T",
       &mm,
       &nn,
       &kk,
       &alpha,
       this->values.data(),
       &kk,
       B.values.data(),
       &nn,
       &beta,
       C.values.data(),
       &mm);
}


template <typename number>
void
LAPACKFullMatrix<number>::TmTmult(FullMatrix<number> &            C,
                                  const LAPACKFullMatrix<number> &B,
                                  const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(B.state));
  Assert(this->m() == B.n(), ExcDimensionMismatch(this->m(), B.n()));
  Assert(C.n() == B.m(), ExcDimensionMismatch(C.n(), B.m()));
  Assert(C.m() == this->n(), ExcDimensionMismatch(this->n(), C.m()));
  const types::blas_int mm    = this->n();
  const types::blas_int nn    = B.m();
  const types::blas_int kk    = B.n();
  const number          alpha = 1.;
  const number          beta  = (adding ? 1. : 0.);

  // since FullMatrix stores the matrix in transposed order compared to this
  // matrix, compute B * A = (A^T * B^T)^T
  gemm("N",
       "N",
       &nn,
       &mm,
       &kk,
       &alpha,
       B.values.data(),
       &nn,
       this->values.data(),
       &kk,
       &beta,
       &C(0, 0),
       &nn);
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_lu_factorization()
{
  Assert(state == matrix, ExcState(state));
  state = LAPACKSupport::unusable;

  const types::blas_int mm     = this->m();
  const types::blas_int nn     = this->n();
  number *const         values = this->values.data();
  ipiv.resize(mm);
  types::blas_int info = 0;
  getrf(&mm, &nn, values, &mm, ipiv.data(), &info);

  Assert(info >= 0, ExcInternalError());

  // if info >= 0, the factorization has been completed
  state = lu;

  AssertThrow(info == 0, LACExceptions::ExcSingular());
}



template <typename number>
void
LAPACKFullMatrix<number>::set_property(const Property p)
{
  property = p;
}



template <typename number>
number
LAPACKFullMatrix<number>::l1_norm() const
{
  const char type('O');
  return norm(type);
}



template <typename number>
number
LAPACKFullMatrix<number>::linfty_norm() const
{
  const char type('I');
  return norm(type);
}



template <typename number>
number
LAPACKFullMatrix<number>::frobenius_norm() const
{
  const char type('F');
  return norm(type);
}



template <typename number>
number
LAPACKFullMatrix<number>::norm(const char type) const
{
  std::lock_guard<std::mutex> lock(mutex);

  Assert(state == LAPACKSupport::matrix ||
           state == LAPACKSupport::inverse_matrix,
         ExcMessage("norms can be called in matrix state only."));

  const types::blas_int N      = this->n();
  const types::blas_int M      = this->m();
  const number *const   values = this->values.data();
  if (property == symmetric)
    {
      const types::blas_int lda = std::max<types::blas_int>(1, N);
      const types::blas_int lwork =
        (type == 'I' || type == 'O') ? std::max<types::blas_int>(1, N) : 0;
      work.resize(lwork);
      return lansy(&type, &LAPACKSupport::L, &N, values, &lda, work.data());
    }
  else
    {
      const types::blas_int lda = std::max<types::blas_int>(1, M);
      const types::blas_int lwork =
        (type == 'I') ? std::max<types::blas_int>(1, M) : 0;
      work.resize(lwork);
      return lange(&type, &M, &N, values, &lda, work.data());
    }
}



template <typename number>
number
LAPACKFullMatrix<number>::trace() const
{
  Assert(state == LAPACKSupport::matrix ||
           state == LAPACKSupport::inverse_matrix,
         ExcMessage("Trace can be called in matrix state only."));
  Assert(this->n() == this->m(), ExcDimensionMismatch(this->n(), this->m()));

  number tr = 0;
  for (size_type i = 0; i < this->m(); ++i)
    tr += (*this)(i, i);

  return tr;
}



template <typename number>
void
LAPACKFullMatrix<number>::compute_cholesky_factorization()
{
  Assert(state == matrix, ExcState(state));
  Assert(property == symmetric, ExcProperty(property));
  state = LAPACKSupport::unusable;

  const types::blas_int mm = this->m();
  const types::blas_int nn = this->n();
  (void)mm;
  Assert(mm == nn, ExcDimensionMismatch(mm, nn));

  number *const         values = this->values.data();
  types::blas_int       info   = 0;
  const types::blas_int lda    = std::max<types::blas_int>(1, nn);
  potrf(&LAPACKSupport::L, &nn, values, &lda, &info);

  // info < 0 : the info-th argument had an illegal value
  Assert(info >= 0, ExcInternalError());

  state = cholesky;
  AssertThrow(info == 0, LACExceptions::ExcSingular());
}



template <typename number>
number
LAPACKFullMatrix<number>::reciprocal_condition_number(const number a_norm) const
{
  std::lock_guard<std::mutex> lock(mutex);
  Assert(state == cholesky, ExcState(state));
  number rcond = 0.;

  const types::blas_int N      = this->m();
  const number *        values = this->values.data();
  types::blas_int       info   = 0;
  const types::blas_int lda    = std::max<types::blas_int>(1, N);
  work.resize(3 * N);
  iwork.resize(N);

  // use the same uplo as in Cholesky
  pocon(&LAPACKSupport::L,
        &N,
        values,
        &lda,
        &a_norm,
        &rcond,
        work.data(),
        iwork.data(),
        &info);

  Assert(info >= 0, ExcInternalError());

  return rcond;
}



template <typename number>
number
LAPACKFullMatrix<number>::reciprocal_condition_number() const
{
  std::lock_guard<std::mutex> lock(mutex);
  Assert(property == upper_triangular || property == lower_triangular,
         ExcProperty(property));
  number rcond = 0.;

  const types::blas_int N      = this->m();
  const number *const   values = this->values.data();
  types::blas_int       info   = 0;
  const types::blas_int lda    = std::max<types::blas_int>(1, N);
  work.resize(3 * N);
  iwork.resize(N);

  const char norm = '1';
  const char diag = 'N';
  const char uplo =
    (property == upper_triangular ? LAPACKSupport::U : LAPACKSupport::L);
  trcon(&norm,
        &uplo,
        &diag,
        &N,
        values,
        &lda,
        &rcond,
        work.data(),
        iwork.data(),
        &info);

  Assert(info >= 0, ExcInternalError());

  return rcond;
}



template <typename number>
void
LAPACKFullMatrix<number>::compute_svd()
{
  Assert(state == matrix, ExcState(state));
  state = LAPACKSupport::unusable;

  const types::blas_int mm = this->m();
  const types::blas_int nn = this->n();
  wr.resize(std::max(mm, nn));
  std::fill(wr.begin(), wr.end(), 0.);
  ipiv.resize(8 * mm);

  svd_u  = std_cxx14::make_unique<LAPACKFullMatrix<number>>(mm, mm);
  svd_vt = std_cxx14::make_unique<LAPACKFullMatrix<number>>(nn, nn);
  types::blas_int info = 0;

  // First determine optimal workspace size
  work.resize(1);
  types::blas_int lwork = -1;

  // TODO double check size
  std::vector<typename numbers::NumberTraits<number>::real_type> real_work;
  if (numbers::NumberTraits<number>::is_complex)
    {
      // This array is only used by the complex versions.
      std::size_t min = std::min(this->m(), this->n());
      std::size_t max = std::max(this->m(), this->n());
      real_work.resize(
        std::max(5 * min * min + 5 * min, 2 * max * min + 2 * min * min + min));
    }

  // make sure that the first entry in the work array is clear, in case the
  // routine does not completely overwrite the memory:
  work[0] = number();
  internal::LAPACKFullMatrixImplementation::gesdd_helper(LAPACKSupport::A,
                                                         mm,
                                                         nn,
                                                         this->values,
                                                         wr,
                                                         svd_u->values,
                                                         svd_vt->values,
                                                         work,
                                                         real_work,
                                                         ipiv,
                                                         lwork,
                                                         info);

  AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("gesdd", info));
  // Resize the work array. Add one to the size computed by LAPACK to be on
  // the safe side.
  lwork = static_cast<types::blas_int>(std::abs(work[0]) + 1);

  work.resize(lwork);
  // Do the actual SVD.
  internal::LAPACKFullMatrixImplementation::gesdd_helper(LAPACKSupport::A,
                                                         mm,
                                                         nn,
                                                         this->values,
                                                         wr,
                                                         svd_u->values,
                                                         svd_vt->values,
                                                         work,
                                                         real_work,
                                                         ipiv,
                                                         lwork,
                                                         info);
  AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("gesdd", info));

  work.resize(0);
  ipiv.resize(0);

  state = LAPACKSupport::svd;
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_inverse_svd(const double threshold)
{
  if (state == LAPACKSupport::matrix)
    compute_svd();

  Assert(state == LAPACKSupport::svd, ExcState(state));

  const typename numbers::NumberTraits<number>::real_type one(1.0);
  const double lim = std::abs(wr[0]) * threshold;
  for (size_type i = 0; i < wr.size(); ++i)
    {
      if (std::abs(wr[i]) > lim)
        wr[i] = one / wr[i];
      else
        wr[i] = 0.;
    }
  state = LAPACKSupport::inverse_svd;
}



template <typename number>
void
LAPACKFullMatrix<number>::compute_inverse_svd_with_kernel(
  const unsigned int kernel_size)
{
  if (state == LAPACKSupport::matrix)
    compute_svd();

  Assert(state == LAPACKSupport::svd, ExcState(state));

  const typename numbers::NumberTraits<number>::real_type one(1.0);
  const unsigned int                                      n_wr = wr.size();
  for (size_type i = 0; i < n_wr - kernel_size; ++i)
    wr[i] = one / wr[i];
  for (size_type i = n_wr - kernel_size; i < n_wr; ++i)
    wr[i] = 0.;
  state = LAPACKSupport::inverse_svd;
}



template <typename number>
void
LAPACKFullMatrix<number>::invert()
{
  Assert(state == matrix || state == lu || state == cholesky, ExcState(state));
  const types::blas_int mm = this->m();
  const types::blas_int nn = this->n();
  Assert(nn == mm, ExcNotQuadratic());

  number *const   values = this->values.data();
  types::blas_int info   = 0;

  if (property != symmetric)
    {
      if (state == matrix)
        compute_lu_factorization();

      ipiv.resize(mm);
      inv_work.resize(mm);
      getri(&mm, values, &mm, ipiv.data(), inv_work.data(), &mm, &info);
    }
  else
    {
      if (state == matrix)
        compute_cholesky_factorization();

      const types::blas_int lda = std::max<types::blas_int>(1, nn);
      potri(&LAPACKSupport::L, &nn, values, &lda, &info);
      // inverse is stored in lower diagonal, set the upper diagonal
      // appropriately:
      for (types::blas_int i = 0; i < nn; ++i)
        for (types::blas_int j = i + 1; j < nn; ++j)
          this->el(i, j) = this->el(j, i);
    }

  Assert(info >= 0, ExcInternalError());
  AssertThrow(info == 0, LACExceptions::ExcSingular());

  state = inverse_matrix;
}



template <typename number>
void
LAPACKFullMatrix<number>::solve(Vector<number> &v, const bool transposed) const
{
  Assert(this->m() == this->n(), LACExceptions::ExcNotQuadratic());
  AssertDimension(this->m(), v.size());
  const char *          trans  = transposed ? &T : &N;
  const types::blas_int nn     = this->n();
  const number *const   values = this->values.data();
  const types::blas_int n_rhs  = 1;
  types::blas_int       info   = 0;

  if (state == lu)
    {
      getrs(
        trans, &nn, &n_rhs, values, &nn, ipiv.data(), v.begin(), &nn, &info);
    }
  else if (state == cholesky)
    {
      potrs(&LAPACKSupport::L, &nn, &n_rhs, values, &nn, v.begin(), &nn, &info);
    }
  else if (property == upper_triangular || property == lower_triangular)
    {
      const char uplo =
        (property == upper_triangular ? LAPACKSupport::U : LAPACKSupport::L);

      const types::blas_int lda = nn;
      const types::blas_int ldb = nn;
      trtrs(
        &uplo, trans, "N", &nn, &n_rhs, values, &lda, v.begin(), &ldb, &info);
    }
  else
    {
      Assert(false,
             ExcMessage(
               "The matrix has to be either factorized or triangular."));
    }

  Assert(info == 0, ExcInternalError());
}



template <typename number>
void
LAPACKFullMatrix<number>::solve(LAPACKFullMatrix<number> &B,
                                const bool                transposed) const
{
  Assert(B.state == matrix, ExcState(B.state));

  Assert(this->m() == this->n(), LACExceptions::ExcNotQuadratic());
  AssertDimension(this->m(), B.m());
  const char *          trans  = transposed ? &T : &N;
  const types::blas_int nn     = this->n();
  const number *const   values = this->values.data();
  const types::blas_int n_rhs  = B.n();
  types::blas_int       info   = 0;

  if (state == lu)
    {
      getrs(trans,
            &nn,
            &n_rhs,
            values,
            &nn,
            ipiv.data(),
            B.values.data(),
            &nn,
            &info);
    }
  else if (state == cholesky)
    {
      potrs(&LAPACKSupport::L,
            &nn,
            &n_rhs,
            values,
            &nn,
            B.values.data(),
            &nn,
            &info);
    }
  else if (property == upper_triangular || property == lower_triangular)
    {
      const char uplo =
        (property == upper_triangular ? LAPACKSupport::U : LAPACKSupport::L);

      const types::blas_int lda = nn;
      const types::blas_int ldb = nn;
      trtrs(&uplo,
            trans,
            "N",
            &nn,
            &n_rhs,
            values,
            &lda,
            B.values.data(),
            &ldb,
            &info);
    }
  else
    {
      Assert(false,
             ExcMessage(
               "The matrix has to be either factorized or triangular."));
    }

  Assert(info == 0, ExcInternalError());
}



template <typename number>
number
LAPACKFullMatrix<number>::determinant() const
{
  Assert(this->m() == this->n(), LACExceptions::ExcNotQuadratic());

  // LAPACK doesn't offer a function to compute a matrix determinant.
  // This is due to the difficulty in maintaining numerical accuracy, as the
  // calculations are likely to overflow or underflow. See
  // http://www.netlib.org/lapack/faq.html#_are_there_routines_in_lapack_to_compute_determinants
  //
  // However, after a PLU decomposition one can compute this by multiplication
  // of the diagonal entries with one another. One must take into consideration
  // the number of permutations (row swaps) stored in the P matrix.
  //
  // See the implementations in the blaze library (detNxN)
  // https://bitbucket.org/blaze-lib/blaze
  // and also
  // https://dualm.wordpress.com/2012/01/06/computing-determinant-in-fortran/
  // http://icl.cs.utk.edu/lapack-forum/viewtopic.php?p=341&#p336
  // for further information.
  Assert(state == lu, ExcState(state));
  Assert(ipiv.size() == this->m(), ExcInternalError());
  number det = 1.0;
  for (size_type i = 0; i < this->m(); ++i)
    det *=
      (ipiv[i] == types::blas_int(i + 1) ? this->el(i, i) : -this->el(i, i));
  return det;
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_eigenvalues(const bool right, const bool left)
{
  Assert(state == matrix, ExcState(state));
  const types::blas_int nn = this->n();
  wr.resize(nn);
  wi.resize(nn);
  if (right)
    vr.resize(nn * nn);
  if (left)
    vl.resize(nn * nn);

  types::blas_int info  = 0;
  types::blas_int lwork = 1;
  const char      jobvr = (right) ? V : N;
  const char      jobvl = (left) ? V : N;

  /*
   * The LAPACK routine xGEEV requires a sufficiently large work array; the
   * minimum requirement is
   *
   * work.size >= 4*nn.
   *
   * However, for better performance, a larger work array may be needed. The
   * first call determines the optimal work size and the second does the work.
   */
  lwork = -1;
  work.resize(1);

  std::vector<typename numbers::NumberTraits<number>::real_type> real_work;
  if (numbers::NumberTraits<number>::is_complex)
    // This array is only used by the complex versions.
    real_work.resize(2 * this->m());
  internal::LAPACKFullMatrixImplementation::geev_helper(jobvl,
                                                        jobvr,
                                                        this->values,
                                                        this->m(),
                                                        wr,
                                                        wi,
                                                        vl,
                                                        vr,
                                                        work,
                                                        real_work,
                                                        lwork,
                                                        info);

  // geev returns info=0 on success. Since we only queried the optimal size
  // for work, everything else would not be acceptable.
  Assert(info == 0, ExcInternalError());
  // Allocate working array according to suggestion (same strategy as was
  // noted in compute_svd).
  lwork = static_cast<types::blas_int>(std::abs(work[0]) + 1);

  // resize workspace array
  work.resize(lwork);

  // Finally compute the eigenvalues.
  internal::LAPACKFullMatrixImplementation::geev_helper(jobvl,
                                                        jobvr,
                                                        this->values,
                                                        this->m(),
                                                        wr,
                                                        wi,
                                                        vl,
                                                        vr,
                                                        work,
                                                        real_work,
                                                        lwork,
                                                        info);

  Assert(info >= 0, ExcInternalError());
  // TODO:[GK] What if the QR method fails?
  if (info != 0)
    std::cerr << "LAPACK error in geev" << std::endl;

  state = LAPACKSupport::State(LAPACKSupport::eigenvalues | unusable);
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_eigenvalues_symmetric(
  const number        lower_bound,
  const number        upper_bound,
  const number        abs_accuracy,
  Vector<number> &    eigenvalues,
  FullMatrix<number> &eigenvectors)
{
  Assert(state == matrix, ExcState(state));
  const types::blas_int nn = (this->n() > 0 ? this->n() : 1);
  Assert(static_cast<size_type>(nn) == this->m(), ExcNotQuadratic());

  wr.resize(nn);
  LAPACKFullMatrix<number> matrix_eigenvectors(nn, nn);

  number *const values_A            = this->values.data();
  number *const values_eigenvectors = matrix_eigenvectors.values.data();

  types::blas_int              info(0), lwork(-1), n_eigenpairs(0);
  const char *const            jobz(&V);
  const char *const            uplo(&U);
  const char *const            range(&V);
  const types::blas_int *const dummy(&one);
  std::vector<types::blas_int> iwork(static_cast<size_type>(5 * nn));
  std::vector<types::blas_int> ifail(static_cast<size_type>(nn));


  /*
   * The LAPACK routine xSYEVX requires a sufficiently large work array; the
   * minimum requirement is
   *
   * work.size >= 8*nn.
   *
   * However, for better performance, a larger work array may be needed. The
   * first call determines the optimal work size and the second does the work.
   */
  work.resize(1);

  syevx(jobz,
        range,
        uplo,
        &nn,
        values_A,
        &nn,
        &lower_bound,
        &upper_bound,
        dummy,
        dummy,
        &abs_accuracy,
        &n_eigenpairs,
        wr.data(),
        values_eigenvectors,
        &nn,
        work.data(),
        &lwork,
        iwork.data(),
        ifail.data(),
        &info);
  // syevx returns info=0 on success. Since we only queried the optimal size
  // for work, everything else would not be acceptable.
  Assert(info == 0, ExcInternalError());
  // Allocate working array according to suggestion (same strategy as was noted
  // in compute_svd).
  lwork = static_cast<types::blas_int>(std::abs(work[0]) + 1);
  work.resize(static_cast<size_type>(lwork));

  // Finally compute the eigenvalues.
  syevx(jobz,
        range,
        uplo,
        &nn,
        values_A,
        &nn,
        &lower_bound,
        &upper_bound,
        dummy,
        dummy,
        &abs_accuracy,
        &n_eigenpairs,
        wr.data(),
        values_eigenvectors,
        &nn,
        work.data(),
        &lwork,
        iwork.data(),
        ifail.data(),
        &info);

  // Negative return value implies a wrong argument. This should be internal.
  Assert(info >= 0, ExcInternalError());
  if (info != 0)
    std::cerr << "LAPACK error in syevx" << std::endl;

  eigenvalues.reinit(n_eigenpairs);
  eigenvectors.reinit(nn, n_eigenpairs, true);

  for (size_type i = 0; i < static_cast<size_type>(n_eigenpairs); ++i)
    {
      eigenvalues(i) = wr[i];
      size_type col_begin(i * nn);
      for (size_type j = 0; j < static_cast<size_type>(nn); ++j)
        {
          eigenvectors(j, i) = values_eigenvectors[col_begin + j];
        }
    }

  state = LAPACKSupport::State(unusable);
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_generalized_eigenvalues_symmetric(
  LAPACKFullMatrix<number> &   B,
  const number                 lower_bound,
  const number                 upper_bound,
  const number                 abs_accuracy,
  Vector<number> &             eigenvalues,
  std::vector<Vector<number>> &eigenvectors,
  const types::blas_int        itype)
{
  Assert(state == matrix, ExcState(state));
  const types::blas_int nn = (this->n() > 0 ? this->n() : 1);
  Assert(static_cast<size_type>(nn) == this->m(), ExcNotQuadratic());
  Assert(B.m() == B.n(), ExcNotQuadratic());
  Assert(static_cast<size_type>(nn) == B.n(), ExcDimensionMismatch(nn, B.n()));

  wr.resize(nn);
  LAPACKFullMatrix<number> matrix_eigenvectors(nn, nn);

  number *const values_A            = this->values.data();
  number *const values_B            = B.values.data();
  number *const values_eigenvectors = matrix_eigenvectors.values.data();

  types::blas_int              info(0), lwork(-1), n_eigenpairs(0);
  const char *const            jobz(&V);
  const char *const            uplo(&U);
  const char *const            range(&V);
  const types::blas_int *const dummy(&one);
  iwork.resize(static_cast<size_type>(5 * nn));
  std::vector<types::blas_int> ifail(static_cast<size_type>(nn));


  /*
   * The LAPACK routine xSYGVX requires a sufficiently large work array; the
   * minimum requirement is
   *
   * work.size >= 8*nn.
   *
   * However, for better performance, a larger work array may be needed. The
   * first call determines the optimal work size and the second does the work.
   */
  work.resize(1);

  sygvx(&itype,
        jobz,
        range,
        uplo,
        &nn,
        values_A,
        &nn,
        values_B,
        &nn,
        &lower_bound,
        &upper_bound,
        dummy,
        dummy,
        &abs_accuracy,
        &n_eigenpairs,
        wr.data(),
        values_eigenvectors,
        &nn,
        work.data(),
        &lwork,
        iwork.data(),
        ifail.data(),
        &info);
  // sygvx returns info=0 on success. Since we only queried the optimal size
  // for work, everything else would not be acceptable.
  Assert(info == 0, ExcInternalError());
  // Allocate working array according to suggestion (same strategy as was
  // noted in compute_svd).
  lwork = static_cast<types::blas_int>(std::abs(work[0]) + 1);

  // resize workspace arrays
  work.resize(static_cast<size_type>(lwork));

  // Finally compute the generalized eigenvalues.
  sygvx(&itype,
        jobz,
        range,
        uplo,
        &nn,
        values_A,
        &nn,
        values_B,
        &nn,
        &lower_bound,
        &upper_bound,
        dummy,
        dummy,
        &abs_accuracy,
        &n_eigenpairs,
        wr.data(),
        values_eigenvectors,
        &nn,
        work.data(),
        &lwork,
        iwork.data(),
        ifail.data(),
        &info);

  // Negative return value implies a wrong argument. This should be internal.
  Assert(info >= 0, ExcInternalError());
  if (info != 0)
    std::cerr << "LAPACK error in sygvx" << std::endl;

  eigenvalues.reinit(n_eigenpairs);
  eigenvectors.resize(n_eigenpairs);

  for (size_type i = 0; i < static_cast<size_type>(n_eigenpairs); ++i)
    {
      eigenvalues(i) = wr[i];
      size_type col_begin(i * nn);
      eigenvectors[i].reinit(nn, true);
      for (size_type j = 0; j < static_cast<size_type>(nn); ++j)
        {
          eigenvectors[i](j) = values_eigenvectors[col_begin + j];
        }
    }

  state = LAPACKSupport::State(unusable);
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_generalized_eigenvalues_symmetric(
  LAPACKFullMatrix<number> &   B,
  std::vector<Vector<number>> &eigenvectors,
  const types::blas_int        itype)
{
  Assert(state == matrix, ExcState(state));
  const types::blas_int nn = this->n();
  Assert(static_cast<size_type>(nn) == this->m(), ExcNotQuadratic());
  Assert(B.m() == B.n(), ExcNotQuadratic());
  Assert(static_cast<size_type>(nn) == B.n(), ExcDimensionMismatch(nn, B.n()));
  Assert(eigenvectors.size() <= static_cast<size_type>(nn),
         ExcMessage("eigenvectors.size() > matrix.n()"));

  wr.resize(nn);
  wi.resize(nn); // This is set purely for consistency reasons with the
  // eigenvalues() function.

  number *const values_A = this->values.data();
  number *const values_B = B.values.data();

  types::blas_int   info  = 0;
  types::blas_int   lwork = -1;
  const char *const jobz  = (eigenvectors.size() > 0) ? (&V) : (&N);
  const char *const uplo  = (&U);

  /*
   * The LAPACK routine xSYGV requires a sufficiently large work array; the
   * minimum requirement is
   *
   * work.size >= 3*nn - 1.
   *
   * However, for better performance, a larger work array may be needed. The
   * first call determines the optimal work size and the second does the work.
   */
  work.resize(1);

  sygv(&itype,
       jobz,
       uplo,
       &nn,
       values_A,
       &nn,
       values_B,
       &nn,
       wr.data(),
       work.data(),
       &lwork,
       &info);
  // sygv returns info=0 on success. Since we only queried the optimal size
  // for work, everything else would not be acceptable.
  Assert(info == 0, ExcInternalError());
  // Allocate working array according to suggestion (same strategy as was
  // noted in compute_svd).
  lwork = static_cast<types::blas_int>(std::abs(work[0]) + 1);

  // resize workspace array
  work.resize(static_cast<size_type>(lwork));

  // Finally compute the generalized eigenvalues.
  sygv(&itype,
       jobz,
       uplo,
       &nn,
       values_A,
       &nn,
       values_B,
       &nn,
       wr.data(),
       work.data(),
       &lwork,
       &info);
  // Negative return value implies a wrong argument. This should be internal.

  Assert(info >= 0, ExcInternalError());
  if (info != 0)
    std::cerr << "LAPACK error in sygv" << std::endl;

  for (size_type i = 0; i < eigenvectors.size(); ++i)
    {
      size_type col_begin(i * nn);
      eigenvectors[i].reinit(nn, true);
      for (size_type j = 0; j < static_cast<size_type>(nn); ++j)
        {
          eigenvectors[i](j) = values_A[col_begin + j];
        }
    }
  state = LAPACKSupport::State(LAPACKSupport::eigenvalues | unusable);
}


template <typename number>
void
LAPACKFullMatrix<number>::print_formatted(std::ostream &     out,
                                          const unsigned int precision,
                                          const bool         scientific,
                                          const unsigned int width_,
                                          const char *       zero_string,
                                          const double       denominator,
                                          const double       threshold) const
{
  unsigned int width = width_;

  Assert((!this->empty()) || (this->n() + this->m() == 0), ExcInternalError());
  Assert(state == LAPACKSupport::matrix ||
           state == LAPACKSupport::inverse_matrix ||
           state == LAPACKSupport::cholesky,
         ExcState(state));

  // set output format, but store old
  // state
  std::ios::fmtflags old_flags     = out.flags();
  std::streamsize    old_precision = out.precision(precision);

  if (scientific)
    {
      out.setf(std::ios::scientific, std::ios::floatfield);
      if (!width)
        width = precision + 7;
    }
  else
    {
      out.setf(std::ios::fixed, std::ios::floatfield);
      if (!width)
        width = precision + 2;
    }

  for (size_type i = 0; i < this->m(); ++i)
    {
      // Cholesky is stored in lower triangular, so just output this part:
      const size_type nc = state == LAPACKSupport::cholesky ? i + 1 : this->n();
      for (size_type j = 0; j < nc; ++j)
        // we might have complex numbers, so use abs also to check for nan
        // since there is no isnan on complex numbers
        if (std::isnan(std::abs((*this)(i, j))))
          out << std::setw(width) << (*this)(i, j) << ' ';
        else if (std::abs(this->el(i, j)) > threshold)
          out << std::setw(width) << this->el(i, j) * denominator << ' ';
        else
          out << std::setw(width) << zero_string << ' ';
      out << std::endl;
    }

  AssertThrow(out, ExcIO());
  // reset output format
  out.flags(old_flags);
  out.precision(old_precision);
}


//----------------------------------------------------------------------//

template <typename number>
void
PreconditionLU<number>::initialize(const LAPACKFullMatrix<number> &M)
{
  matrix = &M;
  mem    = nullptr;
}


template <typename number>
void
PreconditionLU<number>::initialize(const LAPACKFullMatrix<number> &M,
                                   VectorMemory<Vector<number>> &  V)
{
  matrix = &M;
  mem    = &V;
}


template <typename number>
void
PreconditionLU<number>::vmult(Vector<number> &      dst,
                              const Vector<number> &src) const
{
  dst = src;
  matrix->solve(dst, false);
}


template <typename number>
void
PreconditionLU<number>::Tvmult(Vector<number> &      dst,
                               const Vector<number> &src) const
{
  dst = src;
  matrix->solve(dst, true);
}


template <typename number>
void
PreconditionLU<number>::vmult(BlockVector<number> &      dst,
                              const BlockVector<number> &src) const
{
  Assert(mem != nullptr, ExcNotInitialized());
  Vector<number> *aux = mem->alloc();
  *aux                = src;
  matrix->solve(*aux, false);
  dst = *aux;
}


template <typename number>
void
PreconditionLU<number>::Tvmult(BlockVector<number> &      dst,
                               const BlockVector<number> &src) const
{
  Assert(mem != nullptr, ExcNotInitialized());
  Vector<number> *aux = mem->alloc();
  *aux                = src;
  matrix->solve(*aux, true);
  dst = *aux;
}



#include "lapack_full_matrix.inst"


DEAL_II_NAMESPACE_CLOSE
