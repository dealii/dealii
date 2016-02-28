// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2015 by the deal.II authors
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


#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>

#include <iostream>
#include <iomanip>

DEAL_II_NAMESPACE_OPEN

using namespace LAPACKSupport;

template <typename number>
LAPACKFullMatrix<number>::LAPACKFullMatrix (const size_type n)
  :
  TransposeTable<number> (n,n),
  state (matrix)
{}


template <typename number>
LAPACKFullMatrix<number>::LAPACKFullMatrix (const size_type m,
                                            const size_type n)
  :
  TransposeTable<number> (m, n),
  state (matrix)
{}


template <typename number>
LAPACKFullMatrix<number>::LAPACKFullMatrix (const LAPACKFullMatrix &M)
  :
  TransposeTable<number> (M),
  state (matrix)
{}


template <typename number>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator = (const LAPACKFullMatrix<number> &M)
{
  TransposeTable<number>::operator=(M);
  state = LAPACKSupport::matrix;
  return *this;
}


template <typename number>
void
LAPACKFullMatrix<number>::reinit (const size_type n)
{
  this->TransposeTable<number>::reinit (n, n);
  state = LAPACKSupport::matrix;
}


template <typename number>
void
LAPACKFullMatrix<number>::reinit (const size_type m,
                                  const size_type n)
{
  this->TransposeTable<number>::reinit (m, n);
  state = LAPACKSupport::matrix;
}


template <typename number>
template <typename number2>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator = (const FullMatrix<number2> &M)
{
  Assert (this->n_rows() == M.n_rows(), ExcDimensionMismatch(this->n_rows(), M.n_rows()));
  Assert (this->n_cols() == M.n(), ExcDimensionMismatch(this->n_cols(), M.n()));
  for (size_type i=0; i<this->n_rows(); ++i)
    for (size_type j=0; j<this->n_cols(); ++j)
      (*this)(i,j) = M(i,j);

  state = LAPACKSupport::matrix;
  return *this;
}


template <typename number>
template <typename number2>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator = (const SparseMatrix<number2> &M)
{
  Assert (this->n_rows() == M.n(), ExcDimensionMismatch(this->n_rows(), M.n()));
  Assert (this->n_cols() == M.m(), ExcDimensionMismatch(this->n_cols(), M.m()));
  for (size_type i=0; i<this->n_rows(); ++i)
    for (size_type j=0; j<this->n_cols(); ++j)
      (*this)(i,j) = M.el(i,j);

  state = LAPACKSupport::matrix;
  return *this;
}


template <typename number>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator = (const double d)
{
  (void)d;
  Assert (d==0, ExcScalarAssignmentOnlyForZeroValue());

  if (this->n_elements() != 0)
    this->reset_values();

  state = LAPACKSupport::matrix;
  return *this;
}


template <typename number>
void
LAPACKFullMatrix<number>::vmult (
  Vector<number>       &w,
  const Vector<number> &v,
  const bool            adding) const
{
  const int mm = this->n_rows();
  const int nn = this->n_cols();
  const number alpha = 1.;
  const number beta = (adding ? 1. : 0.);
  const number null = 0.;

  switch (state)
    {
    case matrix:
    case inverse_matrix:
    {
      AssertDimension(v.size(), this->n_cols());
      AssertDimension(w.size(), this->n_rows());

      gemv("N", &mm, &nn, &alpha, &this->values[0], &mm, v.val, &one, &beta, w.val, &one);
      break;
    }
    case svd:
    {
      AssertDimension(v.size(), this->n_cols());
      AssertDimension(w.size(), this->n_rows());
      // Compute V^T v
      work.resize(std::max(mm,nn));
      gemv("N", &nn, &nn, &alpha, &svd_vt->values[0], &nn, v.val, &one, &null, &work[0], &one);
      // Multiply by singular values
      for (size_type i=0; i<wr.size(); ++i)
        work[i] *= wr[i];
      // Multiply with U
      gemv("N", &mm, &mm, &alpha, &svd_u->values[0], &mm, &work[0], &one, &beta, w.val, &one);
      break;
    }
    case inverse_svd:
    {
      AssertDimension(w.size(), this->n_cols());
      AssertDimension(v.size(), this->n_rows());
      // Compute U^T v
      work.resize(std::max(mm,nn));
      gemv("T", &mm, &mm, &alpha, &svd_u->values[0], &mm, v.val, &one, &null, &work[0], &one);
      // Multiply by singular values
      for (size_type i=0; i<wr.size(); ++i)
        work[i] *= wr[i];
      // Multiply with V
      gemv("T", &nn, &nn, &alpha, &svd_vt->values[0], &nn, &work[0], &one, &beta, w.val, &one);
      break;
    }
    default:
      Assert (false, ExcState(state));
    }
}


template <typename number>
void
LAPACKFullMatrix<number>::Tvmult (
  Vector<number>       &w,
  const Vector<number> &v,
  const bool            adding) const
{
  const int mm = this->n_rows();
  const int nn = this->n_cols();
  const number alpha = 1.;
  const number beta = (adding ? 1. : 0.);
  const number null = 0.;

  switch (state)
    {
    case matrix:
    case inverse_matrix:
    {
      AssertDimension(w.size(), this->n_cols());
      AssertDimension(v.size(), this->n_rows());

      gemv("T", &mm, &nn, &alpha, &this->values[0], &mm, v.val, &one, &beta, w.val, &one);
      break;
    }
    case svd:
    {
      AssertDimension(w.size(), this->n_cols());
      AssertDimension(v.size(), this->n_rows());

      // Compute U^T v
      work.resize(std::max(mm,nn));
      gemv("T", &mm, &mm, &alpha, &svd_u->values[0], &mm, v.val, &one, &null, &work[0], &one);
      // Multiply by singular values
      for (size_type i=0; i<wr.size(); ++i)
        work[i] *= wr[i];
      // Multiply with V
      gemv("T", &nn, &nn, &alpha, &svd_vt->values[0], &nn, &work[0], &one, &beta, w.val, &one);
      break;
      case inverse_svd:
      {
        AssertDimension(v.size(), this->n_cols());
        AssertDimension(w.size(), this->n_rows());

        // Compute V^T v
        work.resize(std::max(mm,nn));
        gemv("N", &nn, &nn, &alpha, &svd_vt->values[0], &nn, v.val, &one, &null, &work[0], &one);
        // Multiply by singular values
        for (size_type i=0; i<wr.size(); ++i)
          work[i] *= wr[i];
        // Multiply with U
        gemv("N", &mm, &mm, &alpha, &svd_u->values[0], &mm, &work[0], &one, &beta, w.val, &one);
        break;
      }
    }
    default:
      Assert (false, ExcState(state));
    }
}


template <typename number>
void
LAPACKFullMatrix<number>::vmult_add (Vector<number>       &w,
                                     const Vector<number> &v) const
{
  vmult(w, v, true);
}


template <typename number>
void
LAPACKFullMatrix<number>::Tvmult_add (Vector<number>       &w,
                                      const Vector<number> &v) const
{
  Tvmult(w, v, true);
}


template <typename number>
void
LAPACKFullMatrix<number>::mmult(LAPACKFullMatrix<number>       &C,
                                const LAPACKFullMatrix<number> &B,
                                const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(state));
  Assert(C.state == matrix || C.state == inverse_matrix, ExcState(state));
  Assert (this->n_cols() == B.n_rows(), ExcDimensionMismatch(this->n_cols(), B.n_rows()));
  Assert (C.n_cols() == B.n_cols(), ExcDimensionMismatch(C.n_cols(), B.n_cols()));
  Assert (C.n_rows() == this->n_rows(), ExcDimensionMismatch(this->n_rows(), C.n_rows()));
  const int mm = this->n_rows();
  const int nn = B.n_cols();
  const int kk = this->n_cols();
  const number alpha = 1.;
  const number beta = (adding ? 1. : 0.);

  gemm("N", "N", &mm, &nn, &kk, &alpha, &this->values[0], &mm, &B.values[0],
       &kk, &beta, &C.values[0], &mm);
}


template <typename number>
void
LAPACKFullMatrix<number>::mmult(FullMatrix<number>             &C,
                                const LAPACKFullMatrix<number> &B,
                                const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(state));
  Assert (this->n_cols() == B.n_rows(), ExcDimensionMismatch(this->n_cols(), B.n_rows()));
  Assert (C.n_cols() == B.n_cols(), ExcDimensionMismatch(C.n_cols(), B.n_cols()));
  Assert (C.n_rows() == this->n_rows(), ExcDimensionMismatch(this->n_rows(), C.n_rows()));
  const int mm = this->n_rows();
  const int nn = B.n_cols();
  const int kk = this->n_cols();
  const number alpha = 1.;
  const number beta = (adding ? 1. : 0.);

  // since FullMatrix stores the matrix in transposed order compared to this
  // matrix, compute B^T * A^T = (A * B)^T
  gemm("T", "T", &nn, &mm, &kk, &alpha, &B.values[0], &kk, &this->values[0],
       &mm, &beta, &C(0,0), &nn);
}



template <typename number>
void
LAPACKFullMatrix<number>::Tmmult(LAPACKFullMatrix<number>       &C,
                                 const LAPACKFullMatrix<number> &B,
                                 const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(state));
  Assert(C.state == matrix || C.state == inverse_matrix, ExcState(state));
  Assert (this->n_rows() == B.n_rows(), ExcDimensionMismatch(this->n_rows(), B.n_rows()));
  Assert (C.n_cols() == B.n_cols(), ExcDimensionMismatch(C.n_cols(), B.n_cols()));
  Assert (C.n_rows() == this->n_cols(), ExcDimensionMismatch(this->n_cols(), C.n_rows()));
  const int mm = this->n_cols();
  const int nn = B.n_cols();
  const int kk = B.n_rows();
  const number alpha = 1.;
  const number beta = (adding ? 1. : 0.);

  gemm("T", "N", &mm, &nn, &kk, &alpha, &this->values[0], &kk, &B.values[0],
       &kk, &beta, &C.values[0], &mm);
}


template <typename number>
void
LAPACKFullMatrix<number>::Tmmult(FullMatrix<number>             &C,
                                 const LAPACKFullMatrix<number> &B,
                                 const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(state));
  Assert (this->n_rows() == B.n_rows(), ExcDimensionMismatch(this->n_rows(), B.n_rows()));
  Assert (C.n_cols() == B.n_cols(), ExcDimensionMismatch(C.n_cols(), B.n_cols()));
  Assert (C.n_rows() == this->n_cols(), ExcDimensionMismatch(this->n_cols(), C.n_rows()));
  const int mm = this->n_cols();
  const int nn = B.n_cols();
  const int kk = B.n_rows();
  const number alpha = 1.;
  const number beta = (adding ? 1. : 0.);

  // since FullMatrix stores the matrix in transposed order compared to this
  // matrix, compute B^T * A = (A^T * B)^T
  gemm("T", "N", &nn, &mm, &kk, &alpha, &B.values[0], &kk, &this->values[0],
       &kk, &beta, &C(0,0), &nn);
}


template <typename number>
void
LAPACKFullMatrix<number>::mTmult(LAPACKFullMatrix<number>       &C,
                                 const LAPACKFullMatrix<number> &B,
                                 const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(state));
  Assert(C.state == matrix || C.state == inverse_matrix, ExcState(state));
  Assert (this->n_cols() == B.n_cols(), ExcDimensionMismatch(this->n_cols(), B.n_cols()));
  Assert (C.n_cols() == B.n_rows(), ExcDimensionMismatch(C.n_cols(), B.n_rows()));
  Assert (C.n_rows() == this->n_rows(), ExcDimensionMismatch(this->n_rows(), C.n_rows()));
  const int mm = this->n_rows();
  const int nn = B.n_rows();
  const int kk = B.n_cols();
  const number alpha = 1.;
  const number beta = (adding ? 1. : 0.);

  gemm("N", "T", &mm, &nn, &kk, &alpha, &this->values[0], &mm, &B.values[0],
       &nn, &beta, &C.values[0], &mm);
}



template <typename number>
void
LAPACKFullMatrix<number>::mTmult(FullMatrix<number>             &C,
                                 const LAPACKFullMatrix<number> &B,
                                 const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(state));
  Assert (this->n_cols() == B.n_cols(), ExcDimensionMismatch(this->n_cols(), B.n_cols()));
  Assert (C.n_cols() == B.n_rows(), ExcDimensionMismatch(C.n_cols(), B.n_rows()));
  Assert (C.n_rows() == this->n_rows(), ExcDimensionMismatch(this->n_rows(), C.n_rows()));
  const int mm = this->n_rows();
  const int nn = B.n_rows();
  const int kk = B.n_cols();
  const number alpha = 1.;
  const number beta = (adding ? 1. : 0.);

  // since FullMatrix stores the matrix in transposed order compared to this
  // matrix, compute B * A^T = (A * B^T)^T
  gemm("N", "T", &nn, &mm, &kk, &alpha, &B.values[0], &nn, &this->values[0],
       &mm, &beta, &C(0,0), &nn);
}


template <typename number>
void
LAPACKFullMatrix<number>::TmTmult(LAPACKFullMatrix<number>       &C,
                                  const LAPACKFullMatrix<number> &B,
                                  const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(state));
  Assert(C.state == matrix || C.state == inverse_matrix, ExcState(state));
  Assert (this->n_rows() == B.n_cols(), ExcDimensionMismatch(this->n_rows(), B.n_cols()));
  Assert (C.n_cols() == B.n_rows(), ExcDimensionMismatch(C.n_cols(), B.n_rows()));
  Assert (C.n_rows() == this->n_cols(), ExcDimensionMismatch(this->n_cols(), C.n_rows()));
  const int mm = this->n_cols();
  const int nn = B.n_rows();
  const int kk = B.n_cols();
  const number alpha = 1.;
  const number beta = (adding ? 1. : 0.);

  gemm("T", "T", &mm, &nn, &kk, &alpha, &this->values[0], &kk, &B.values[0],
       &nn, &beta, &C.values[0], &mm);
}


template <typename number>
void
LAPACKFullMatrix<number>::TmTmult(FullMatrix<number>             &C,
                                  const LAPACKFullMatrix<number> &B,
                                  const bool                      adding) const
{
  Assert(state == matrix || state == inverse_matrix, ExcState(state));
  Assert(B.state == matrix || B.state == inverse_matrix, ExcState(state));
  Assert (this->n_rows() == B.n_cols(), ExcDimensionMismatch(this->n_rows(), B.n_cols()));
  Assert (C.n_cols() == B.n_rows(), ExcDimensionMismatch(C.n_cols(), B.n_rows()));
  Assert (C.n_rows() == this->n_cols(), ExcDimensionMismatch(this->n_cols(), C.n_rows()));
  const int mm = this->n_cols();
  const int nn = B.n_rows();
  const int kk = B.n_cols();
  const number alpha = 1.;
  const number beta = (adding ? 1. : 0.);

  // since FullMatrix stores the matrix in transposed order compared to this
  // matrix, compute B * A = (A^T * B^T)^T
  gemm("N", "N", &nn, &mm, &kk, &alpha, &B.values[0], &nn, &this->values[0],
       &kk, &beta, &C(0,0), &nn);
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_lu_factorization()
{
  Assert(state == matrix, ExcState(state));
  const int mm = this->n_rows();
  const int nn = this->n_cols();
  number *values = const_cast<number *> (&this->values[0]);
  ipiv.resize(mm);
  int info = 0;
  getrf(&mm, &nn, values, &mm, &ipiv[0], &info);

  AssertThrow(info >= 0, ExcInternalError());
  AssertThrow(info == 0, LACExceptions::ExcSingular());

  state = lu;
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_svd()
{
  Assert(state == matrix, ExcState(state));
  state = LAPACKSupport::unusable;

  const int mm = this->n_rows();
  const int nn = this->n_cols();
  number *values = const_cast<number *> (&this->values[0]);
  wr.resize(std::max(mm,nn));
  std::fill(wr.begin(), wr.end(), 0.);
  ipiv.resize(8*mm);

  svd_u.reset (new LAPACKFullMatrix<number>(mm,mm));
  svd_vt.reset (new LAPACKFullMatrix<number>(nn,nn));
  number *mu  = const_cast<number *> (&svd_u->values[0]);
  number *mvt = const_cast<number *> (&svd_vt->values[0]);
  int info = 0;

  // First determine optimal workspace size
  work.resize(1);
  int lwork = -1;
  gesdd(&LAPACKSupport::A, &mm, &nn, values, &mm,
        &wr[0], mu, &mm, mvt, &nn,
        &work[0], &lwork, &ipiv[0], &info);
  AssertThrow (info==0, LAPACKSupport::ExcErrorCode("gesdd", info));
  // Resize the work array. Add one to the size computed by LAPACK to be on
  // the safe side.
  lwork = static_cast<int>(work[0] + 1);

  work.resize(lwork);
  // Do the actual SVD.
  gesdd(&LAPACKSupport::A, &mm, &nn, values, &mm,
        &wr[0], mu, &mm, mvt, &nn,
        &work[0], &lwork, &ipiv[0], &info);
  AssertThrow (info==0, LAPACKSupport::ExcErrorCode("gesdd", info));

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

  Assert (state==LAPACKSupport::svd, ExcState(state));

  const double lim = wr[0]*threshold;
  for (size_type i=0; i<wr.size(); ++i)
    {
      if (wr[i] > lim)
        wr[i] = 1./wr[i];
      else
        wr[i] = 0.;
    }
  state = LAPACKSupport::inverse_svd;
}


template <typename number>
void
LAPACKFullMatrix<number>::invert()
{
  Assert(state == matrix || state == lu,
         ExcState(state));
  const int mm = this->n_rows();
  const int nn = this->n_cols();
  Assert (nn == mm, ExcNotQuadratic());

  number *values = const_cast<number *> (&this->values[0]);
  ipiv.resize(mm);
  int info = 0;

  if (state == matrix)
    {
      getrf(&mm, &nn, values, &mm, &ipiv[0], &info);

      AssertThrow(info >= 0, ExcInternalError());
      AssertThrow(info == 0, LACExceptions::ExcSingular());
    }

  inv_work.resize (mm);
  getri(&mm, values, &mm, &ipiv[0], &inv_work[0], &mm, &info);

  AssertThrow(info >= 0, ExcInternalError());
  AssertThrow(info == 0, LACExceptions::ExcSingular());

  state = inverse_matrix;
}


template <typename number>
void
LAPACKFullMatrix<number>::apply_lu_factorization(Vector<number> &v,
                                                 const bool transposed) const
{
  Assert(state == lu, ExcState(state));
  Assert(this->n_rows() == this->n_cols(),
         LACExceptions::ExcNotQuadratic());
  AssertDimension(this->n_rows(), v.size());

  const char *trans = transposed ? &T : &N;
  const int nn = this->n_cols();
  const number *values = &this->values[0];
  int info = 0;

  getrs(trans, &nn, &one, values, &nn, &ipiv[0],
        v.begin(), &nn, &info);

  AssertThrow(info == 0, ExcInternalError());
}


template <typename number>
void
LAPACKFullMatrix<number>::apply_lu_factorization(LAPACKFullMatrix<number> &B,
                                                 const bool transposed) const
{
  Assert(state == lu, ExcState(state));
  Assert(B.state == matrix, ExcState(state));
  Assert(this->n_rows() == this->n_cols(), LACExceptions::ExcNotQuadratic());
  AssertDimension(this->n_rows(), B.n_rows());

  const char *trans = transposed ? &T : &N;
  const int nn = this->n_cols();
  const int kk = B.n_cols();
  const number *values = &this->values[0];
  int info = 0;

  getrs(trans, &nn, &kk, values, &nn, &ipiv[0], &B.values[0], &nn, &info);

  AssertThrow(info == 0, ExcInternalError());
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_eigenvalues(const bool right,
                                              const bool left)
{
  Assert(state == matrix, ExcState(state));
  const int nn = this->n_cols();
  wr.resize(nn);
  wi.resize(nn);
  if (right) vr.resize(nn*nn);
  if (left)  vl.resize(nn*nn);

  number *values = const_cast<number *> (&this->values[0]);

  int info  = 0;
  int lwork = 1;
  const char *const jobvr = (right) ? (&V) : (&N);
  const char *const jobvl = (left)  ? (&V) : (&N);

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

  geev(jobvl, jobvr, &nn, values, &nn,
       &wr[0], &wi[0],
       &vl[0], &nn, &vr[0], &nn,
       &work[0], &lwork, &info);
  // geev returns info=0 on success. Since we only queried the optimal size
  // for work, everything else would not be acceptable.
  Assert (info == 0, ExcInternalError());
  // Allocate working array according to suggestion (same strategy as was
  // noted in compute_svd).
  lwork = static_cast<int>(work[0] + 1);

  // resize workspace array
  work.resize((size_type ) lwork);

  // Finally compute the eigenvalues.
  geev(jobvl, jobvr, &nn, values, &nn,
       &wr[0], &wi[0],
       &vl[0], &nn, &vr[0], &nn,
       &work[0], &lwork, &info);
  // Negative return value implies a wrong argument. This should be internal.

  Assert (info >=0, ExcInternalError());
//TODO:[GK] What if the QR method fails?
  if (info != 0)
    std::cerr << "LAPACK error in geev" << std::endl;

  state = LAPACKSupport::State(eigenvalues | unusable);
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_eigenvalues_symmetric(const number        lower_bound,
                                                        const number        upper_bound,
                                                        const number        abs_accuracy,
                                                        Vector<number>     &eigenvalues,
                                                        FullMatrix<number> &eigenvectors)
{
  Assert(state == matrix, ExcState(state));
  const int nn = (this->n_cols() > 0 ? this->n_cols() : 1);
  Assert(static_cast<size_type>(nn) == this->n_rows(), ExcNotQuadratic());

  wr.resize(nn);
  LAPACKFullMatrix<number> matrix_eigenvectors(nn, nn);

  number *values_A = const_cast<number *> (&this->values[0]);
  number *values_eigenvectors = const_cast<number *> (&matrix_eigenvectors.values[0]);

  int info(0),
      lwork(-1),
      n_eigenpairs(0);
  const char *const jobz(&V);
  const char *const uplo(&U);
  const char *const range(&V);
  const int *const  dummy(&one);
  std::vector<int> iwork(static_cast<size_type> (5*nn));
  std::vector<int> ifail(static_cast<size_type> (nn));


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

  syevx (jobz, range,
         uplo, &nn, values_A, &nn,
         &lower_bound, &upper_bound,
         dummy, dummy, &abs_accuracy,
         &n_eigenpairs, &wr[0], values_eigenvectors,
         &nn, &work[0], &lwork, &iwork[0],
         &ifail[0], &info);
  // syevx returns info=0 on success. Since we only queried the optimal size
  // for work, everything else would not be acceptable.
  Assert (info == 0, ExcInternalError());
  // Allocate working array according to suggestion (same strategy as was noted in
  // compute_svd).
  lwork = static_cast<int>(work[0] + 1);
  work.resize(static_cast<size_type> (lwork));

  // Finally compute the eigenvalues.
  syevx (jobz, range,
         uplo, &nn, values_A, &nn,
         &lower_bound, &upper_bound,
         dummy, dummy, &abs_accuracy,
         &n_eigenpairs, &wr[0], values_eigenvectors,
         &nn, &work[0], &lwork, &iwork[0],
         &ifail[0], &info);

  // Negative return value implies a wrong argument. This should be internal.
  Assert (info >=0, ExcInternalError());
  if (info != 0)
    std::cerr << "LAPACK error in syevx" << std::endl;

  eigenvalues.reinit(n_eigenpairs);
  eigenvectors.reinit(nn, n_eigenpairs, true);

  for (size_type i=0; i < static_cast<size_type> (n_eigenpairs); ++i)
    {
      eigenvalues(i) = wr[i];
      size_type col_begin(i*nn);
      for (size_type j=0; j < static_cast<size_type> (nn); ++j)
        {
          eigenvectors(j,i) = values_eigenvectors[col_begin+j];
        }
    }

  state = LAPACKSupport::State(unusable);
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_generalized_eigenvalues_symmetric(
  LAPACKFullMatrix<number> &B,
  const number lower_bound,
  const number upper_bound,
  const number abs_accuracy,
  Vector<number> &eigenvalues,
  std::vector<Vector<number> > &eigenvectors,
  const int itype)
{
  Assert(state == matrix, ExcState(state));
  const int nn = (this->n_cols() > 0 ? this->n_cols() : 1);
  Assert(static_cast<size_type>(nn) == this->n_rows(), ExcNotQuadratic());
  Assert(B.n_rows() == B.n_cols(), ExcNotQuadratic());
  Assert(static_cast<size_type>(nn) == B.n_cols(),
         ExcDimensionMismatch (nn, B.n_cols()));

  wr.resize(nn);
  LAPACKFullMatrix<number> matrix_eigenvectors(nn, nn);

  number *values_A = const_cast<number *> (&this->values[0]);
  number *values_B = const_cast<number *> (&B.values[0]);
  number *values_eigenvectors = const_cast<number *> (&matrix_eigenvectors.values[0]);

  int info(0),
      lwork(-1),
      n_eigenpairs(0);
  const char *const jobz(&V);
  const char *const uplo(&U);
  const char *const range(&V);
  const int *const  dummy(&one);
  std::vector<int> iwork(static_cast<size_type> (5*nn));
  std::vector<int> ifail(static_cast<size_type> (nn));


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

  sygvx (&itype, jobz, range, uplo, &nn, values_A, &nn,
         values_B, &nn, &lower_bound, &upper_bound,
         dummy, dummy, &abs_accuracy, &n_eigenpairs,
         &wr[0], values_eigenvectors, &nn, &work[0],
         &lwork, &iwork[0], &ifail[0], &info);
  // sygvx returns info=0 on success. Since we only queried the optimal size
  // for work, everything else would not be acceptable.
  Assert (info == 0, ExcInternalError());
  // Allocate working array according to suggestion (same strategy as was
  // noted in compute_svd).
  lwork = static_cast<int>(work[0] + 1);

  // resize workspace arrays
  work.resize(static_cast<size_type> (lwork));

  // Finally compute the generalized eigenvalues.
  sygvx (&itype, jobz, range, uplo, &nn, values_A, &nn,
         values_B, &nn, &lower_bound, &upper_bound,
         dummy, dummy, &abs_accuracy, &n_eigenpairs,
         &wr[0], values_eigenvectors, &nn, &work[0],
         &lwork, &iwork[0], &ifail[0], &info);

  // Negative return value implies a wrong argument. This should be internal.
  Assert (info >=0, ExcInternalError());
  if (info != 0)
    std::cerr << "LAPACK error in sygvx" << std::endl;

  eigenvalues.reinit(n_eigenpairs);
  eigenvectors.resize(n_eigenpairs);

  for (size_type i=0; i < static_cast<size_type> (n_eigenpairs); ++i)
    {
      eigenvalues(i) = wr[i];
      size_type col_begin(i*nn);
      eigenvectors[i].reinit(nn, true);
      for (size_type j=0; j < static_cast<size_type> (nn); ++j)
        {
          eigenvectors[i](j) = values_eigenvectors[col_begin+j];
        }
    }

  state = LAPACKSupport::State(unusable);
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_generalized_eigenvalues_symmetric (
  LAPACKFullMatrix<number> &B,
  std::vector<Vector<number> > &eigenvectors,
  const int itype)
{
  Assert(state == matrix, ExcState(state));
  const int nn = this->n_cols();
  Assert(static_cast<size_type>(nn) == this->n_rows(), ExcNotQuadratic());
  Assert(B.n_rows() == B.n_cols(), ExcNotQuadratic());
  Assert(static_cast<size_type>(nn) == B.n_cols(),
         ExcDimensionMismatch (nn, B.n_cols()));
  Assert(eigenvectors.size() <= static_cast<size_type>(nn),
         ExcMessage ("eigenvectors.size() > matrix.n_cols()"));

  wr.resize(nn);
  wi.resize(nn); //This is set purely for consistency reasons with the
  //eigenvalues() function.

  number *values_A = const_cast<number *> (&this->values[0]);
  number *values_B = const_cast<number *> (&B.values[0]);

  int info  = 0;
  int lwork = -1;
  const char *const jobz = (eigenvectors.size() > 0) ? (&V) : (&N);
  const char *const uplo = (&U);

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

  sygv (&itype, jobz, uplo, &nn, values_A, &nn,
        values_B, &nn,
        &wr[0], &work[0], &lwork, &info);
  // sygv returns info=0 on success. Since we only queried the optimal size
  // for work, everything else would not be acceptable.
  Assert (info == 0, ExcInternalError());
  // Allocate working array according to suggestion (same strategy as was
  // noted in compute_svd).
  lwork = static_cast<int>(work[0] + 1);

  // resize workspace array
  work.resize(static_cast<size_type>(lwork));

  // Finally compute the generalized eigenvalues.
  sygv (&itype, jobz, uplo, &nn, values_A, &nn,
        values_B, &nn,
        &wr[0], &work[0], &lwork, &info);
  // Negative return value implies a wrong argument. This should be internal.

  Assert (info >=0, ExcInternalError());
  if (info != 0)
    std::cerr << "LAPACK error in sygv" << std::endl;

  for (size_type i=0; i < eigenvectors.size(); ++i)
    {
      size_type col_begin(i*nn);
      eigenvectors[i].reinit(nn, true);
      for (size_type j=0; j < static_cast<size_type>(nn); ++j)
        {
          eigenvectors[i](j) = values_A[col_begin+j];
        }
    }
  state = LAPACKSupport::State(eigenvalues | unusable);
}


template <typename number>
void
LAPACKFullMatrix<number>::print_formatted (
  std::ostream       &out,
  const unsigned int  precision,
  const bool          scientific,
  const unsigned int  width_,
  const char         *zero_string,
  const double        denominator,
  const double        threshold) const
{
  unsigned int width = width_;

  Assert ((!this->empty()) || (this->n_cols()+this->n_rows()==0),
          ExcInternalError());

  // set output format, but store old
  // state
  std::ios::fmtflags old_flags = out.flags();
  unsigned int old_precision = out.precision (precision);

  if (scientific)
    {
      out.setf (std::ios::scientific, std::ios::floatfield);
      if (!width)
        width = precision+7;
    }
  else
    {
      out.setf (std::ios::fixed, std::ios::floatfield);
      if (!width)
        width = precision+2;
    }

  for (size_type i=0; i<this->n_rows(); ++i)
    {
      for (size_type j=0; j<this->n_cols(); ++j)
        if (std::fabs(this->el(i,j)) > threshold)
          out << std::setw(width)
              << this->el(i,j) * denominator << ' ';
        else
          out << std::setw(width) << zero_string << ' ';
      out << std::endl;
    };

  AssertThrow (out, ExcIO());
  // reset output format
  out.flags (old_flags);
  out.precision(old_precision);
}


//----------------------------------------------------------------------//

template <typename number>
void
PreconditionLU<number>::initialize(const LAPACKFullMatrix<number> &M)
{
  matrix = &M;
  mem = 0;
}


template <typename number>
void
PreconditionLU<number>::initialize(const LAPACKFullMatrix<number> &M,
                                   VectorMemory<Vector<number> > &V)
{
  matrix = &M;
  mem = &V;
}


template <typename number>
void
PreconditionLU<number>::vmult(Vector<number> &dst,
                              const Vector<number> &src) const
{
  dst = src;
  matrix->apply_lu_factorization(dst, false);
}


template <typename number>
void
PreconditionLU<number>::Tvmult(Vector<number> &dst,
                               const Vector<number> &src) const
{
  dst = src;
  matrix->apply_lu_factorization(dst, true);
}


template <typename number>
void
PreconditionLU<number>::vmult(BlockVector<number> &dst,
                              const BlockVector<number> &src) const
{
  Assert(mem != 0, ExcNotInitialized());
  Vector<number> *aux = mem->alloc();
  *aux = src;
  matrix->apply_lu_factorization(*aux, false);
  dst = *aux;
}


template <typename number>
void
PreconditionLU<number>::Tvmult(BlockVector<number> &dst,
                               const BlockVector<number> &src) const
{
  Assert(mem != 0, ExcNotInitialized());
  Vector<number> *aux = mem->alloc();
  *aux = src;
  matrix->apply_lu_factorization(*aux, true);
  dst = *aux;
}



#include "lapack_full_matrix.inst"


DEAL_II_NAMESPACE_CLOSE
