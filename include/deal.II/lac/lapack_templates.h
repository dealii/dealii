// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


#ifndef __LAPACK_TEMPLATES_H
#define __LAPACK_TEMPLATES_H

#include <deal.II/base/config.h>
#include <deal.II/lac/lapack_support.h>

extern "C"
{
// vector update of the form y += alpha*x with a scalar, x,y vectors
  void daxpy_ (const int *n, const double *alpha, const double *x,
               const int *incx, double *y, const int *incy);
  void saxpy_ (const int *n, const float *alpha, const float *x,
               const int *incx, float *y, const int *incy);
// General Matrix
// Matrix vector product
  void dgemv_ (const char *trans, const int *m, const int *n,
               const double *alpha, const double *A, const int *lda,
               const double *x, const int *incx,
               const double *b, double *y, const int *incy);
  void sgemv_ (const char *trans, const int *m, const int *n,
               const float *alpha, const float *A, const int *lda,
               const float *x, const int *incx,
               const float *b, float *y, const int *incy);
// Matrix matrix product
  void dgemm_ (const char *transa, const char *transb,
               const int *m, const int *n, const int *k,
               const double *alpha, const double *A, const int *lda,
               const double *B, const int *ldb,
               const double *beta, double *C, const int *ldc);
  void sgemm_ (const char *transa, const char *transb,
               const int *m, const int *n, const int *k,
               const float *alpha, const float *A, const int *lda,
               const float *B, const int *ldb,
               const float *beta, float *C, const int *ldc);
// Compute LU factorization
  void dgetrf_ (const int *m, const int *n, double *A,
                const int *lda, int *ipiv, int *info);
  void sgetrf_ (const int *m, const int *n, float *A,
                const int *lda, int *ipiv, int *info);
// Apply forward/backward substitution to LU factorization
  void dgetrs_ (const char *trans, const int *n, const int *nrhs,
                const double *A, const int *lda, const int *ipiv,
                double *b, const int *ldb, int *info);
  void sgetrs_ (const char *trans, const int *n, const int *nrhs,
                const float *A, const int *lda, const int *ipiv,
                float *b, const int *ldb, int *info);
// Invert matrix from LU factorization
  void dgetri_ (const int *n, double *A, const int *lda,
                int *ipiv, double *inv_work, const int *lwork, int *info);
  void sgetri_ (const int *n, float *A, const int *lda,
                int *ipiv, float *inv_work, const int *lwork, int *info);
// Compute QR factorization (Householder)
  void dgeqrf_ (const int *m, const int *n, double *A,
                const int *lda, double *tau, double *work,
                const int *lwork, int *info);
  void sgeqrf_ (const int *m, const int *n, float *A,
                const int *lda, float *tau, float *work,
                const int *lwork, int *info);
// Compute vector Q^T B, where Q is the result from dgeqrf_
  void dormqr_ (const char *side, const char *trans, const int *m,
                const int *n, const int *k, const double *A, const int *lda,
                const double *tau, double *B, const int *ldb,
                double *work, const int *lwork, int *info);
  void sormqr_ (const char *side, const char *trans, const int *m,
                const int *n, const int *k, const float *A, const int *lda,
                const float *tau, float *B, const int *ldb,
                float *work, const int *lwork, int *info);
// Compute matrix Q from the result of dgeqrf_
  void dorgqr_ (const int *m, const int *n, const int *k, const double *A,
                const int *lda, const double *tau, double *work, const int *lwork,
                int *info);
  void sorgqr_ (const int *m, const int *n, const int *k, const float *A,
                const int *lda, const float *tau, float *work, const int *lwork,
                int *info);
// Compute Rx = b
  void dtrtrs_ (const char *uplo, const char *trans,
                const char *diag, const int *n, const int *n_rhs,
                const double *A, const int *lda, double *B, const int *ldb,
                int *info);
  void strtrs_ (const char *uplo, const char *trans,
                const char *diag, const int *n, const int *n_rhs,
                const float *A, const int *lda, float *B, const int *ldb,
                int *info);
// Compute eigenvalues and vectors
  void dgeev_ (const char *jobvl, const char *jobvr,
               const int *n, double *A, const int *lda,
               double *lambda_re, double *lambda_im,
               double *vl, const int *ldvl,
               double *vr, const int *ldva,
               double *work, const int *lwork,
               int *info);
  void sgeev_ (const char *jobvl, const char *jobvr,
               const int *n, float *A, const int *lda,
               float *lambda_re, float *lambda_im,
               float *vl, const int *ldvl,
               float *vr, const int *ldva,
               float *work, const int *lwork,
               int *info);
// Compute eigenvalues and vectors (expert)
  void dgeevx_ (const char *balanc, const char *jobvl, const char *jobvr,
                const char *sense,
                const int *n, double *A, const int *lda,
                double *lambda_re, double *lambda_im,
                double *vl, const int *ldvl,
                double *vr, const int *ldvr,
                int *ilo, int *ihi,
                double *scale, double *abnrm,
                double *rconde, double *rcondv,
                double *work, const int *lwork,
                int *iwork, int *info);
  void sgeevx_ (const char *balanc, const char *jobvl, const char *jobvr,
                const char *sense,
                const int *n, float *A, const int *lda,
                float *lambda_re, float *lambda_im,
                float *vl, const int *ldvl,
                float *vr, const int *ldvr,
                int *ilo, int *ihi,
                float *scale, float *abnrm,
                float *rconde, float *rcondv,
                float *work, const int *lwork,
                int *iwork, int *info);
// Eigenvalues for a symmetric matrix
  void dsyev_ (const char *jobz, const char *uplo, const int *n,
               double *A, const int *lda, double *w,
               double *work, const int *lwork, int *info);
  void ssyev_ (const char *jobz, const char *uplo, const int *n,
               float *A, const int *lda, float *w,
               float *work, const int *lwork, int *info);
// Same functionality as dsyev_ but with more options: E.g.
// Compute only eigenvalues in a specific interval,
// Compute only eigenvalues with a specific index,
// Set tolerance for eigenvalue computation
  void dsyevx_ (const char *jobz, const char *range,
                const char *uplo, const int *n, double *A, const int *lda,
                const double *vl, const double *vu,
                const int *il, const int *iu, const double *abstol,
                int *m, double *w, double *z,
                const int *ldz, double *work, const int *lwork, int *iwork,
                int *ifail, int *info);
  void ssyevx_ (const char *jobz, const char *range,
                const char *uplo, const int *n, float *A, const int *lda,
                const float *vl, const float *vu,
                const int *il, const int *iu, const float *abstol,
                int *m, float *w, float *z,
                const int *ldz, float *work, const int *lwork, int *iwork,
                int *ifail, int *info);
// Generalized eigenvalues and eigenvectors of
// 1: A*x = lambda*B*x; 2: A*B*x = lambda*x; 3: B*A*x = lambda*x
// A and B are symmetric and B is definite
  void dsygv_ (const int *itype, const char *jobz, const char *uplo,
               const int *n, double *A, const int *lda, double *B,
               const int *ldb, double *w, double *work,
               const int *lwork, int *info);
  void ssygv_ (const int *itype, const char *jobz, const char *uplo,
               const int *n, float *A, const int *lda, float *B,
               const int *ldb, float *w, float *work,
               const int *lwork, int *info);
// Same functionality as dsygv_ but with more options: E.g.
// Compute only eigenvalues in a specific interval,
// Compute only eigenvalues with a specific index,
// Set tolerance for eigenvalue computation
  void dsygvx_ (const int *itype, const char *jobz, const char *range,
                const char *uplo, const int *n, double *A, const int *lda,
                double *B, const int *ldb, const double *vl, const double *vu,
                const int *il, const int *iu, const double *abstol,
                int *m, double *w, double *z,
                const int *ldz, double *work, const int *lwork, int *iwork,
                int *ifail, int *info);
  void ssygvx_ (const int *itype, const char *jobz, const char *range,
                const char *uplo, const int *n, float *A, const int *lda,
                float *B, const int *ldb, const float *vl, const float *vu,
                const int *il, const int *iu, const float *abstol,
                int *m, float *w, float *z,
                const int *ldz, float *work, const int *lwork, int *iwork,
                int *ifail, int *info);
// Compute singular value decomposition using divide and conquer
  void dgesdd_ (const char *jobz,
                const int *m, const int *n, double *A, const int *lda,
                double *s,
                double *u, const int *ldu,
                double *vt, const int *ldvt,
                double *work, const int *lwork,
                int *iwork,
                int *info);
  void sgesdd_ (const char *jobz,
                const int *m, const int *n, float *A, const int *lda,
                float *s,
                float *u, const int *ldu,
                float *vt, const int *ldvt,
                float *work, const int *lwork,
                int *iwork,
                int *info);
// Compute singular value decomposition
  void dgesvd_ (int *jobu, int *jobvt,
                const int *n, const int *m, double *A, const int *lda,
                double *s,
                double *u, const int *ldu,
                double *vt, const int *ldvt,
                double *work, const int *lwork,
                int *info);
  void sgesvd_ (int *jobu, int *jobvt,
                const int *n, const int *m, float *A, const int *lda,
                float *s,
                float *u, const int *ldu,
                float *vt, const int *ldvt,
                float *work, const int *lwork,
                int *info);
// Solve a least squares problem using SVD
  void dgelsd_ (const int *m, const int *n, const int *nrhs,
                const double *A, const int *lda,
                double *B, const int *ldb,
                double *s, const double *rcond,
                int *rank,
                double *work, const int *lwork, int *iwork,
                int *info);
  void sgelsd_ (const int *m, const int *n, const int *nrhs,
                const float *A, const int *lda,
                float *B, const int *ldb,
                float *s, const float *rcond,
                int *rank,
                float *work, const int *lwork, int *iwork,
                int *info);
// Symmetric tridiagonal matrix
  void dstev_ (const char *jobz, const int *n,
               double *d, double *e, double *z,
               const int *ldz, double *work,
               int *info);
  void sstev_ (const char *jobz, const int *n,
               float *d, float *e, float *z,
               const int *ldz, float *work,
               int *info);

}

DEAL_II_NAMESPACE_OPEN


/// Template wrapper for LAPACK functions daxpy and saxpy
template<typename number1, typename number2, typename number3>
inline void
axpy (const int *, const number1 *, const number2 *, const int *, number3 *, const int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DAXPY_
inline void
axpy (const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy)
{
  daxpy_ (n,alpha,x,incx,y,incy);
}
#else
inline void
axpy (const int *, const double *, const double *, const int *, double *, const int *)
{
  Assert (false, LAPACKSupport::ExcMissing("daxpy"));
}
#endif


#ifdef HAVE_SAXPY_
inline void
axpy (const int *n, const float *alpha, const float *x, const int *incx, float *y, const int *incy)
{
  saxpy_ (n,alpha,x,incx,y,incy);
}
#else
inline void
axpy (const int *, const float *, const float *, const int *, float *, const int *)
{
  Assert (false, LAPACKSupport::ExcMissing("saxpy"));
}
#endif


/// Template wrapper for LAPACK functions dgemv and sgemv
template<typename number1, typename number2, typename number3, typename number4, typename number5>
inline void
gemv (const char *, const int *, const int *, const number1 *, const number2 *, const int *, const number3 *, const int *, const number4 *, number5 *, const int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DGEMV_
inline void
gemv (const char *trans, const int *m, const int *n, const double *alpha, const double *A, const int *lda, const double *x, const int *incx, const double *b, double *y, const int *incy)
{
  dgemv_ (trans,m,n,alpha,A,lda,x,incx,b,y,incy);
}
#else
inline void
gemv (const char *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dgemv"));
}
#endif


#ifdef HAVE_SGEMV_
inline void
gemv (const char *trans, const int *m, const int *n, const float *alpha, const float *A, const int *lda, const float *x, const int *incx, const float *b, float *y, const int *incy)
{
  sgemv_ (trans,m,n,alpha,A,lda,x,incx,b,y,incy);
}
#else
inline void
gemv (const char *, const int *, const int *, const float *, const float *, const int *, const float *, const int *, const float *, float *, const int *)
{
  Assert (false, LAPACKSupport::ExcMissing("sgemv"));
}
#endif


/// Template wrapper for LAPACK functions dgemm and sgemm
template<typename number1, typename number2, typename number3, typename number4, typename number5>
inline void
gemm (const char *, const char *, const int *, const int *, const int *, const number1 *, const number2 *, const int *, const number3 *, const int *, const number4 *, number5 *, const int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DGEMM_
inline void
gemm (const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha, const double *A, const int *lda, const double *B, const int *ldb, const double *beta, double *C, const int *ldc)
{
  dgemm_ (transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
}
#else
inline void
gemm (const char *, const char *, const int *, const int *, const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dgemm"));
}
#endif


#ifdef HAVE_SGEMM_
inline void
gemm (const char *transa, const char *transb, const int *m, const int *n, const int *k, const float *alpha, const float *A, const int *lda, const float *B, const int *ldb, const float *beta, float *C, const int *ldc)
{
  sgemm_ (transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
}
#else
inline void
gemm (const char *, const char *, const int *, const int *, const int *, const float *, const float *, const int *, const float *, const int *, const float *, float *, const int *)
{
  Assert (false, LAPACKSupport::ExcMissing("sgemm"));
}
#endif


/// Template wrapper for LAPACK functions dgetrf and sgetrf
template<typename number1>
inline void
getrf (const int *, const int *, number1 *, const int *, int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DGETRF_
inline void
getrf (const int *m, const int *n, double *A, const int *lda, int *ipiv, int *info)
{
  dgetrf_ (m,n,A,lda,ipiv,info);
}
#else
inline void
getrf (const int *, const int *, double *, const int *, int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dgetrf"));
}
#endif


#ifdef HAVE_SGETRF_
inline void
getrf (const int *m, const int *n, float *A, const int *lda, int *ipiv, int *info)
{
  sgetrf_ (m,n,A,lda,ipiv,info);
}
#else
inline void
getrf (const int *, const int *, float *, const int *, int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("sgetrf"));
}
#endif


/// Template wrapper for LAPACK functions dgetrs and sgetrs
template<typename number1, typename number2>
inline void
getrs (const char *, const int *, const int *, const number1 *, const int *, const int *, number2 *, const int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DGETRS_
inline void
getrs (const char *trans, const int *n, const int *nrhs, const double *A, const int *lda, const int *ipiv, double *b, const int *ldb, int *info)
{
  dgetrs_ (trans,n,nrhs,A,lda,ipiv,b,ldb,info);
}
#else
inline void
getrs (const char *, const int *, const int *, const double *, const int *, const int *, double *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dgetrs"));
}
#endif


#ifdef HAVE_SGETRS_
inline void
getrs (const char *trans, const int *n, const int *nrhs, const float *A, const int *lda, const int *ipiv, float *b, const int *ldb, int *info)
{
  sgetrs_ (trans,n,nrhs,A,lda,ipiv,b,ldb,info);
}
#else
inline void
getrs (const char *, const int *, const int *, const float *, const int *, const int *, float *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("sgetrs"));
}
#endif


/// Template wrapper for LAPACK functions dgetri and sgetri
template<typename number1, typename number2>
inline void
getri (const int *, number1 *, const int *, int *, number2 *, const int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DGETRI_
inline void
getri (const int *n, double *A, const int *lda, int *ipiv, double *inv_work, const int *lwork, int *info)
{
  dgetri_ (n,A,lda,ipiv,inv_work,lwork,info);
}
#else
inline void
getri (const int *, double *, const int *, int *, double *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dgetri"));
}
#endif


#ifdef HAVE_SGETRI_
inline void
getri (const int *n, float *A, const int *lda, int *ipiv, float *inv_work, const int *lwork, int *info)
{
  sgetri_ (n,A,lda,ipiv,inv_work,lwork,info);
}
#else
inline void
getri (const int *, float *, const int *, int *, float *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("sgetri"));
}
#endif


/// Template wrapper for LAPACK functions dgeqrf and sgeqrf
template<typename number1, typename number2, typename number3>
inline void
geqrf (const int *, const int *, number1 *, const int *, number2 *, number3 *, const int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DGEQRF_
inline void
geqrf (const int *m, const int *n, double *A, const int *lda, double *tau, double *work, const int *lwork, int *info)
{
  dgeqrf_ (m,n,A,lda,tau,work,lwork,info);
}
#else
inline void
geqrf (const int *, const int *, double *, const int *, double *, double *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dgeqrf"));
}
#endif


#ifdef HAVE_SGEQRF_
inline void
geqrf (const int *m, const int *n, float *A, const int *lda, float *tau, float *work, const int *lwork, int *info)
{
  sgeqrf_ (m,n,A,lda,tau,work,lwork,info);
}
#else
inline void
geqrf (const int *, const int *, float *, const int *, float *, float *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("sgeqrf"));
}
#endif


/// Template wrapper for LAPACK functions dormqr and sormqr
template<typename number1, typename number2, typename number3, typename number4>
inline void
ormqr (const char *, const char *, const int *, const int *, const int *, const number1 *, const int *, const number2 *, number3 *, const int *, number4 *, const int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DORMQR_
inline void
ormqr (const char *side, const char *trans, const int *m, const int *n, const int *k, const double *A, const int *lda, const double *tau, double *B, const int *ldb, double *work, const int *lwork, int *info)
{
  dormqr_ (side,trans,m,n,k,A,lda,tau,B,ldb,work,lwork,info);
}
#else
inline void
ormqr (const char *, const char *, const int *, const int *, const int *, const double *, const int *, const double *, double *, const int *, double *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dormqr"));
}
#endif


#ifdef HAVE_SORMQR_
inline void
ormqr (const char *side, const char *trans, const int *m, const int *n, const int *k, const float *A, const int *lda, const float *tau, float *B, const int *ldb, float *work, const int *lwork, int *info)
{
  sormqr_ (side,trans,m,n,k,A,lda,tau,B,ldb,work,lwork,info);
}
#else
inline void
ormqr (const char *, const char *, const int *, const int *, const int *, const float *, const int *, const float *, float *, const int *, float *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("sormqr"));
}
#endif


/// Template wrapper for LAPACK functions dorgqr and sorgqr
template<typename number1, typename number2, typename number3>
inline void
orgqr (const int *, const int *, const int *, const number1 *, const int *, const number2 *, number3 *, const int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DORGQR_
inline void
orgqr (const int *m, const int *n, const int *k, const double *A, const int *lda, const double *tau, double *work, const int *lwork, int *info)
{
  dorgqr_ (m,n,k,A,lda,tau,work,lwork,info);
}
#else
inline void
orgqr (const int *, const int *, const int *, const double *, const int *, const double *, double *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dorgqr"));
}
#endif


#ifdef HAVE_SORGQR_
inline void
orgqr (const int *m, const int *n, const int *k, const float *A, const int *lda, const float *tau, float *work, const int *lwork, int *info)
{
  sorgqr_ (m,n,k,A,lda,tau,work,lwork,info);
}
#else
inline void
orgqr (const int *, const int *, const int *, const float *, const int *, const float *, float *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("sorgqr"));
}
#endif


/// Template wrapper for LAPACK functions dtrtrs and strtrs
template<typename number1, typename number2>
inline void
trtrs (const char *, const char *, const char *, const int *, const int *, const number1 *, const int *, number2 *, const int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DTRTRS_
inline void
trtrs (const char *uplo, const char *trans, const char *diag, const int *n, const int *n_rhs, const double *A, const int *lda, double *B, const int *ldb, int *info)
{
  dtrtrs_ (uplo,trans,diag,n,n_rhs,A,lda,B,ldb,info);
}
#else
inline void
trtrs (const char *, const char *, const char *, const int *, const int *, const double *, const int *, double *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dtrtrs"));
}
#endif


#ifdef HAVE_STRTRS_
inline void
trtrs (const char *uplo, const char *trans, const char *diag, const int *n, const int *n_rhs, const float *A, const int *lda, float *B, const int *ldb, int *info)
{
  strtrs_ (uplo,trans,diag,n,n_rhs,A,lda,B,ldb,info);
}
#else
inline void
trtrs (const char *, const char *, const char *, const int *, const int *, const float *, const int *, float *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("strtrs"));
}
#endif


/// Template wrapper for LAPACK functions dgeev and sgeev
template<typename number1, typename number2, typename number3, typename number4, typename number5, typename number6>
inline void
geev (const char *, const char *, const int *, number1 *, const int *, number2 *, number3 *, number4 *, const int *, number5 *, const int *, number6 *, const int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DGEEV_
inline void
geev (const char *jobvl, const char *jobvr, const int *n, double *A, const int *lda, double *lambda_re, double *lambda_im, double *vl, const int *ldvl, double *vr, const int *ldva, double *work, const int *lwork, int *info)
{
  dgeev_ (jobvl,jobvr,n,A,lda,lambda_re,lambda_im,vl,ldvl,vr,ldva,work,lwork,info);
}
#else
inline void
geev (const char *, const char *, const int *, double *, const int *, double *, double *, double *, const int *, double *, const int *, double *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dgeev"));
}
#endif


#ifdef HAVE_SGEEV_
inline void
geev (const char *jobvl, const char *jobvr, const int *n, float *A, const int *lda, float *lambda_re, float *lambda_im, float *vl, const int *ldvl, float *vr, const int *ldva, float *work, const int *lwork, int *info)
{
  sgeev_ (jobvl,jobvr,n,A,lda,lambda_re,lambda_im,vl,ldvl,vr,ldva,work,lwork,info);
}
#else
inline void
geev (const char *, const char *, const int *, float *, const int *, float *, float *, float *, const int *, float *, const int *, float *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("sgeev"));
}
#endif


/// Template wrapper for LAPACK functions dgeevx and sgeevx
template<typename number1, typename number2, typename number3, typename number4, typename number5, typename number6, typename number7, typename number8, typename number9, typename number10>
inline void
geevx (const char *, const char *, const char *, const char *, const int *, number1 *, const int *, number2 *, number3 *, number4 *, const int *, number5 *, const int *, int *, int *, number6 *, number7 *, number8 *, number9 *, number10 *, const int *, int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DGEEVX_
inline void
geevx (const char *balanc, const char *jobvl, const char *jobvr, const char *sense, const int *n, double *A, const int *lda, double *lambda_re, double *lambda_im, double *vl, const int *ldvl, double *vr, const int *ldvr, int *ilo, int *ihi, double *scale, double *abnrm, double *rconde, double *rcondv, double *work, const int *lwork, int *iwork, int *info)
{
  dgeevx_ (balanc,jobvl,jobvr,sense,n,A,lda,lambda_re,lambda_im,vl,ldvl,vr,ldvr,ilo,ihi,scale,abnrm,rconde,rcondv,work,lwork,iwork,info);
}
#else
inline void
geevx (const char *, const char *, const char *, const char *, const int *, double *, const int *, double *, double *, double *, const int *, double *, const int *, int *, int *, double *, double *, double *, double *, double *, const int *, int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dgeevx"));
}
#endif


#ifdef HAVE_SGEEVX_
inline void
geevx (const char *balanc, const char *jobvl, const char *jobvr, const char *sense, const int *n, float *A, const int *lda, float *lambda_re, float *lambda_im, float *vl, const int *ldvl, float *vr, const int *ldvr, int *ilo, int *ihi, float *scale, float *abnrm, float *rconde, float *rcondv, float *work, const int *lwork, int *iwork, int *info)
{
  sgeevx_ (balanc,jobvl,jobvr,sense,n,A,lda,lambda_re,lambda_im,vl,ldvl,vr,ldvr,ilo,ihi,scale,abnrm,rconde,rcondv,work,lwork,iwork,info);
}
#else
inline void
geevx (const char *, const char *, const char *, const char *, const int *, float *, const int *, float *, float *, float *, const int *, float *, const int *, int *, int *, float *, float *, float *, float *, float *, const int *, int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("sgeevx"));
}
#endif


/// Template wrapper for LAPACK functions dsyev and ssyev
template<typename number1, typename number2, typename number3>
inline void
syev (const char *, const char *, const int *, number1 *, const int *, number2 *, number3 *, const int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DSYEV_
inline void
syev (const char *jobz, const char *uplo, const int *n, double *A, const int *lda, double *w, double *work, const int *lwork, int *info)
{
  dsyev_ (char *jobz,char *uplo,int *n,double*A,int *lda,double*w,double*work,int *lwork,int *info);
}
#else
inline void
syev (const char *, const char *, const int *, double *, const int *, double *, double *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dsyev"));
}
#endif


#ifdef HAVE_SSYEV_
inline void
syev (const char *jobz, const char *uplo, const int *n, float *A, const int *lda, float *w, float *work, const int *lwork, int *info)
{
  ssyev_ (char *jobz,char *uplo,int *n,double*A,int *lda,double*w,double*work,int *lwork,int *info);
}
#else
inline void
syev (const char *, const char *, const int *, float *, const int *, float *, float *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("ssyev"));
}
#endif


/// Template wrapper for LAPACK functions dsyevx and ssyevx
template<typename number1, typename number2, typename number3, typename number4, typename number5, typename number6, typename number7>
inline void
syevx (const char *, const char *, const char *, const int *, number1 *, const int *, const number2 *, const number3 *, const int *, const int *, const number4 *, int *, number5 *, number6 *, const int *, number7 *, const int *, int *, int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DSYEVX_
inline void
syevx (const char *jobz, const char *range, const char *uplo, const int *n, double *A, const int *lda, const double *vl, const double *vu, const int *il, const int *iu, const double *abstol, int *m, double *w, double *z, const int *ldz, double *work, const int *lwork, int *iwork, int *ifail, int *info)
{
  dsyevx_ (jobz,range,uplo,n,A,lda,vl,vu,il,iu,abstol,m,w,z,ldz,work,lwork,iwork,ifail,info);
}
#else
inline void
syevx (const char *, const char *, const char *, const int *, double *, const int *, const double *, const double *, const int *, const int *, const double *, int *, double *, double *, const int *, double *, const int *, int *, int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dsyevx"));
}
#endif


#ifdef HAVE_SSYEVX_
inline void
syevx (const char *jobz, const char *range, const char *uplo, const int *n, float *A, const int *lda, const float *vl, const float *vu, const int *il, const int *iu, const float *abstol, int *m, float *w, float *z, const int *ldz, float *work, const int *lwork, int *iwork, int *ifail, int *info)
{
  ssyevx_ (jobz,range,uplo,n,A,lda,vl,vu,il,iu,abstol,m,w,z,ldz,work,lwork,iwork,ifail,info);
}
#else
inline void
syevx (const char *, const char *, const char *, const int *, float *, const int *, const float *, const float *, const int *, const int *, const float *, int *, float *, float *, const int *, float *, const int *, int *, int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("ssyevx"));
}
#endif


/// Template wrapper for LAPACK functions dsygv and ssygv
template<typename number1, typename number2, typename number3, typename number4>
inline void
sygv (const int *, const char *, const char *, const int *, number1 *, const int *, number2 *, const int *, number3 *, number4 *, const int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DSYGV_
inline void
sygv (const int *itype, const char *jobz, const char *uplo, const int *n, double *A, const int *lda, double *B, const int *ldb, double *w, double *work, const int *lwork, int *info)
{
  dsygv_ (itype,jobz,uplo,n,A,lda,B,ldb,w,work,lwork,info);
}
#else
inline void
sygv (const int *, const char *, const char *, const int *, double *, const int *, double *, const int *, double *, double *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dsygv"));
}
#endif


#ifdef HAVE_SSYGV_
inline void
sygv (const int *itype, const char *jobz, const char *uplo, const int *n, float *A, const int *lda, float *B, const int *ldb, float *w, float *work, const int *lwork, int *info)
{
  ssygv_ (itype,jobz,uplo,n,A,lda,B,ldb,w,work,lwork,info);
}
#else
inline void
sygv (const int *, const char *, const char *, const int *, float *, const int *, float *, const int *, float *, float *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("ssygv"));
}
#endif


/// Template wrapper for LAPACK functions dsygvx and ssygvx
template<typename number1, typename number2, typename number3, typename number4, typename number5, typename number6, typename number7, typename number8>
inline void
sygvx (const int *, const char *, const char *, const char *, const int *, number1 *, const int *, number2 *, const int *, const number3 *, const number4 *, const int *, const int *, const number5 *, int *, number6 *, number7 *, const int *, number8 *, const int *, int *, int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DSYGVX_
inline void
sygvx (const int *itype, const char *jobz, const char *range, const char *uplo, const int *n, double *A, const int *lda, double *B, const int *ldb, const double *vl, const double *vu, const int *il, const int *iu, const double *abstol, int *m, double *w, double *z, const int *ldz, double *work, const int *lwork, int *iwork, int *ifail, int *info)
{
  dsygvx_ (itype,jobz,range,uplo,n,A,lda,B,ldb,vl,vu,il,iu,abstol,m,w,z,ldz,work,lwork,iwork,ifail,info);
}
#else
inline void
sygvx (const int *, const char *, const char *, const char *, const int *, double *, const int *, double *, const int *, const double *, const double *, const int *, const int *, const double *, int *, double *, double *, const int *, double *, const int *, int *, int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dsygvx"));
}
#endif


#ifdef HAVE_SSYGVX_
inline void
sygvx (const int *itype, const char *jobz, const char *range, const char *uplo, const int *n, float *A, const int *lda, float *B, const int *ldb, const float *vl, const float *vu, const int *il, const int *iu, const float *abstol, int *m, float *w, float *z, const int *ldz, float *work, const int *lwork, int *iwork, int *ifail, int *info)
{
  ssygvx_ (itype,jobz,range,uplo,n,A,lda,B,ldb,vl,vu,il,iu,abstol,m,w,z,ldz,work,lwork,iwork,ifail,info);
}
#else
inline void
sygvx (const int *, const char *, const char *, const char *, const int *, float *, const int *, float *, const int *, const float *, const float *, const int *, const int *, const float *, int *, float *, float *, const int *, float *, const int *, int *, int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("ssygvx"));
}
#endif


/// Template wrapper for LAPACK functions dgesdd and sgesdd
template<typename number1, typename number2, typename number3, typename number4, typename number5>
inline void
gesdd (const char *, const int *, const int *, number1 *, const int *, number2 *, number3 *, const int *, number4 *, const int *, number5 *, const int *, int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DGESDD_
inline void
gesdd (const char *jobz, const int *m, const int *n, double *A, const int *lda, double *s, double *u, const int *ldu, double *vt, const int *ldvt, double *work, const int *lwork, int *iwork, int *info)
{
  dgesdd_ (jobz,m,n,A,lda,s,u,ldu,vt,ldvt,work,lwork,iwork,info);
}
#else
inline void
gesdd (const char *, const int *, const int *, double *, const int *, double *, double *, const int *, double *, const int *, double *, const int *, int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dgesdd"));
}
#endif


#ifdef HAVE_SGESDD_
inline void
gesdd (const char *jobz, const int *m, const int *n, float *A, const int *lda, float *s, float *u, const int *ldu, float *vt, const int *ldvt, float *work, const int *lwork, int *iwork, int *info)
{
  sgesdd_ (jobz,m,n,A,lda,s,u,ldu,vt,ldvt,work,lwork,iwork,info);
}
#else
inline void
gesdd (const char *, const int *, const int *, float *, const int *, float *, float *, const int *, float *, const int *, float *, const int *, int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("sgesdd"));
}
#endif


/// Template wrapper for LAPACK functions dgesvd and sgesvd
template<typename number1, typename number2, typename number3, typename number4, typename number5>
inline void
gesvd (int *, int *, const int *, const int *, number1 *, const int *, number2 *, number3 *, const int *, number4 *, const int *, number5 *, const int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DGESVD_
inline void
gesvd (int *jobu, int *jobvt, const int *n, const int *m, double *A, const int *lda, double *s, double *u, const int *ldu, double *vt, const int *ldvt, double *work, const int *lwork, int *info)
{
  dgesvd_ (jobu,jobvt,n,m,A,lda,s,u,ldu,vt,ldvt,work,lwork,info);
}
#else
inline void
gesvd (int *, int *, const int *, const int *, double *, const int *, double *, double *, const int *, double *, const int *, double *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dgesvd"));
}
#endif


#ifdef HAVE_SGESVD_
inline void
gesvd (int *jobu, int *jobvt, const int *n, const int *m, float *A, const int *lda, float *s, float *u, const int *ldu, float *vt, const int *ldvt, float *work, const int *lwork, int *info)
{
  sgesvd_ (jobu,jobvt,n,m,A,lda,s,u,ldu,vt,ldvt,work,lwork,info);
}
#else
inline void
gesvd (int *, int *, const int *, const int *, float *, const int *, float *, float *, const int *, float *, const int *, float *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("sgesvd"));
}
#endif


/// Template wrapper for LAPACK functions dgelsd and sgelsd
template<typename number1, typename number2, typename number3, typename number4, typename number5>
inline void
gelsd (const int *, const int *, const int *, const number1 *, const int *, number2 *, const int *, number3 *, const number4 *, int *, number5 *, const int *, int *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DGELSD_
inline void
gelsd (const int *m, const int *n, const int *nrhs, const double *A, const int *lda, double *B, const int *ldb, double *s, const double *rcond, int *rank, double *work, const int *lwork, int *iwork, int *info)
{
  dgelsd_ (m,n,nrhs,A,lda,B,ldb,s,rcond,rank,work,lwork,iwork,info);
}
#else
inline void
gelsd (const int *, const int *, const int *, const double *, const int *, double *, const int *, double *, const double *, int *, double *, const int *, int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dgelsd"));
}
#endif


#ifdef HAVE_SGELSD_
inline void
gelsd (const int *m, const int *n, const int *nrhs, const float *A, const int *lda, float *B, const int *ldb, float *s, const float *rcond, int *rank, float *work, const int *lwork, int *iwork, int *info)
{
  sgelsd_ (m,n,nrhs,A,lda,B,ldb,s,rcond,rank,work,lwork,iwork,info);
}
#else
inline void
gelsd (const int *, const int *, const int *, const float *, const int *, float *, const int *, float *, const float *, int *, float *, const int *, int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("sgelsd"));
}
#endif


/// Template wrapper for LAPACK functions dstev and sstev
template<typename number1, typename number2, typename number3, typename number4>
inline void
stev (const char *, const int *, number1 *, number2 *, number3 *, const int *, number4 *, int *)
{
  Assert (false, ExcNotImplemented());
}

#ifdef HAVE_DSTEV_
inline void
stev (const char *jobz, const int *n, double *d, double *e, double *z, const int *ldz, double *work, int *info)
{
  dstev_ (jobz,n,d,e,z,ldz,work,info);
}
#else
inline void
stev (const char *, const int *, double *, double *, double *, const int *, double *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("dstev"));
}
#endif


#ifdef HAVE_SSTEV_
inline void
stev (const char *jobz, const int *n, float *d, float *e, float *z, const int *ldz, float *work, int *info)
{
  sstev_ (jobz,n,d,e,z,ldz,work,info);
}
#else
inline void
stev (const char *, const int *, float *, float *, float *, const int *, float *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("sstev"));
}
#endif


DEAL_II_NAMESPACE_CLOSE

#endif
