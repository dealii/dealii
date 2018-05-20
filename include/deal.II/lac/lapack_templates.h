// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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

#ifndef dealii_lapack_templates_h
#define dealii_lapack_templates_h

#include <deal.II/base/config.h>
#include <deal.II/lac/lapack_support.h>

#ifdef DEAL_II_HAVE_FP_EXCEPTIONS
#  include <cfenv>
#endif

extern "C"
{
  // vector update of the form y += alpha*x with a scalar, x,y vectors
  void
  daxpy_(const dealii::types::blas_int* n,
         const double*                  alpha,
         const double*                  x,
         const dealii::types::blas_int* incx,
         double*                        y,
         const dealii::types::blas_int* incy);
  void
  saxpy_(const dealii::types::blas_int* n,
         const float*                   alpha,
         const float*                   x,
         const dealii::types::blas_int* incx,
         float*                         y,
         const dealii::types::blas_int* incy);
  // General Matrix
  // Matrix vector product
  void
  dgemv_(const char*                    trans,
         const dealii::types::blas_int* m,
         const dealii::types::blas_int* n,
         const double*                  alpha,
         const double*                  A,
         const dealii::types::blas_int* lda,
         const double*                  x,
         const dealii::types::blas_int* incx,
         const double*                  b,
         double*                        y,
         const dealii::types::blas_int* incy);
  void
  sgemv_(const char*                    trans,
         const dealii::types::blas_int* m,
         const dealii::types::blas_int* n,
         const float*                   alpha,
         const float*                   A,
         const dealii::types::blas_int* lda,
         const float*                   x,
         const dealii::types::blas_int* incx,
         const float*                   b,
         float*                         y,
         const dealii::types::blas_int* incy);
  void
  dtrmv_(const char*                    uplo,
         const char*                    trans,
         const char*                    diag,
         const dealii::types::blas_int* N,
         const double*                  A,
         const dealii::types::blas_int* lda,
         double*                        x,
         const dealii::types::blas_int* incx);
  void
  strmv_(const char*                    uplo,
         const char*                    trans,
         const char*                    diag,
         const dealii::types::blas_int* N,
         const float*                   A,
         const dealii::types::blas_int* lda,
         float*                         x,
         const dealii::types::blas_int* incx);
  // Matrix matrix product
  void
  dgemm_(const char*                    transa,
         const char*                    transb,
         const dealii::types::blas_int* m,
         const dealii::types::blas_int* n,
         const dealii::types::blas_int* k,
         const double*                  alpha,
         const double*                  A,
         const dealii::types::blas_int* lda,
         const double*                  B,
         const dealii::types::blas_int* ldb,
         const double*                  beta,
         double*                        C,
         const dealii::types::blas_int* ldc);
  void
  sgemm_(const char*                    transa,
         const char*                    transb,
         const dealii::types::blas_int* m,
         const dealii::types::blas_int* n,
         const dealii::types::blas_int* k,
         const float*                   alpha,
         const float*                   A,
         const dealii::types::blas_int* lda,
         const float*                   B,
         const dealii::types::blas_int* ldb,
         const float*                   beta,
         float*                         C,
         const dealii::types::blas_int* ldc);
  // Symmetric rank-k update
  void
  dsyrk_(const char*                    uplo,
         const char*                    trans,
         const dealii::types::blas_int* n,
         const dealii::types::blas_int* k,
         const double*                  alpha,
         const double*                  A,
         const dealii::types::blas_int* lda,
         const double*                  beta,
         double*                        C,
         const dealii::types::blas_int* ldc);
  void
  ssyrk_(const char*                    uplo,
         const char*                    trans,
         const dealii::types::blas_int* n,
         const dealii::types::blas_int* k,
         const float*                   alpha,
         const float*                   A,
         const dealii::types::blas_int* lda,
         const float*                   beta,
         float*                         C,
         const dealii::types::blas_int* ldc);
  // Compute LU factorization
  void
  dgetrf_(const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          double*                        A,
          const dealii::types::blas_int* lda,
          dealii::types::blas_int*       ipiv,
          dealii::types::blas_int*       info);
  void
  sgetrf_(const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          float*                         A,
          const dealii::types::blas_int* lda,
          dealii::types::blas_int*       ipiv,
          dealii::types::blas_int*       info);
  // Apply forward/backward substitution to LU factorization
  void
  dgetrs_(const char*                    trans,
          const dealii::types::blas_int* n,
          const dealii::types::blas_int* nrhs,
          const double*                  A,
          const dealii::types::blas_int* lda,
          const dealii::types::blas_int* ipiv,
          double*                        b,
          const dealii::types::blas_int* ldb,
          dealii::types::blas_int*       info);
  void
  sgetrs_(const char*                    trans,
          const dealii::types::blas_int* n,
          const dealii::types::blas_int* nrhs,
          const float*                   A,
          const dealii::types::blas_int* lda,
          const dealii::types::blas_int* ipiv,
          float*                         b,
          const dealii::types::blas_int* ldb,
          dealii::types::blas_int*       info);
  // Invert matrix from LU factorization
  void
  dgetri_(const dealii::types::blas_int* n,
          double*                        A,
          const dealii::types::blas_int* lda,
          dealii::types::blas_int*       ipiv,
          double*                        inv_work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       info);
  void
  sgetri_(const dealii::types::blas_int* n,
          float*                         A,
          const dealii::types::blas_int* lda,
          dealii::types::blas_int*       ipiv,
          float*                         inv_work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       info);
  // Compute Cholesky factorization of SPD
  void
  dpotrf_(const char*                    uplo,
          const dealii::types::blas_int* n,
          double*                        A,
          const dealii::types::blas_int* lda,
          dealii::types::blas_int*       info);
  void
  spotrf_(const char*                    uplo,
          const dealii::types::blas_int* n,
          float*                         A,
          const dealii::types::blas_int* lda,
          dealii::types::blas_int*       info);
  // Apply forward/backward substitution to Cholesky factorization
  void
  dpotrs_(const char*                    uplo,
          const dealii::types::blas_int* n,
          const dealii::types::blas_int* nrhs,
          const double*                  A,
          const dealii::types::blas_int* lda,
          double*                        B,
          const dealii::types::blas_int* ldb,
          dealii::types::blas_int*       info);
  void
  spotrs_(const char*                    uplo,
          const dealii::types::blas_int* n,
          const dealii::types::blas_int* nrhs,
          const float*                   A,
          const dealii::types::blas_int* lda,
          float*                         B,
          const dealii::types::blas_int* ldb,
          dealii::types::blas_int*       info);
  // Estimate the reciprocal of the condition number in 1-norm from Cholesky
  void
  dpocon_(const char*                    uplo,
          const dealii::types::blas_int* n,
          const double*                  A,
          const dealii::types::blas_int* lda,
          const double*                  anorm,
          double*                        rcond,
          double*                        work,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       info);
  void
  spocon_(const char*                    uplo,
          const dealii::types::blas_int* n,
          const float*                   A,
          const dealii::types::blas_int* lda,
          const float*                   anorm,
          float*                         rcond,
          float*                         work,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       info);
  // Estimate the reciprocal of the condition number of triangular matrices
  // http://www.netlib.org/lapack/explore-html/da/dba/group__double_o_t_h_e_rcomputational_gaff914510b1673e90752c095f5b9dcedf.html#gaff914510b1673e90752c095f5b9dcedf
  void
  dtrcon_(const char*                    norm,
          const char*                    uplo,
          const char*                    diag,
          const dealii::types::blas_int* n,
          const double*                  A,
          const dealii::types::blas_int* lda,
          double*                        rcond,
          double*                        work,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       info);
  void
  strcon_(const char*                    norm,
          const char*                    uplo,
          const char*                    diag,
          const dealii::types::blas_int* n,
          const float*                   A,
          const dealii::types::blas_int* lda,
          float*                         rcond,
          float*                         work,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       info);
  // Computes the inverse from Cholesky
  void
  dpotri_(const char*                    uplo,
          const dealii::types::blas_int* n,
          double*                        A,
          const dealii::types::blas_int* lda,
          dealii::types::blas_int*       info);
  void
  spotri_(const char*                    uplo,
          const dealii::types::blas_int* n,
          float*                         A,
          const dealii::types::blas_int* lda,
          dealii::types::blas_int*       info);
  // Norms
  double
  dlansy_(const char*                    norm,
          const char*                    uplo,
          const dealii::types::blas_int* n,
          const double*                  A,
          const dealii::types::blas_int* lda,
          double*                        work);
  float
  slansy_(const char*                    norm,
          const char*                    uplo,
          const dealii::types::blas_int* n,
          const float*                   A,
          const dealii::types::blas_int* lda,
          float*                         work);
  double
  dlange_(const char*                    norm,
          const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          const double*                  A,
          const dealii::types::blas_int* lda,
          double*                        work);
  float
  slange_(const char*                    norm,
          const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          const float*                   A,
          const dealii::types::blas_int* lda,
          float*                         work);

  // Compute QR factorization (Householder)
  void
  dgeqrf_(const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          double*                        A,
          const dealii::types::blas_int* lda,
          double*                        tau,
          double*                        work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       info);
  void
  sgeqrf_(const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          float*                         A,
          const dealii::types::blas_int* lda,
          float*                         tau,
          float*                         work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       info);
  // Compute vector Q^T B, where Q is the result from dgeqrf_
  void
  dormqr_(const char*                    side,
          const char*                    trans,
          const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          const dealii::types::blas_int* k,
          const double*                  A,
          const dealii::types::blas_int* lda,
          const double*                  tau,
          double*                        B,
          const dealii::types::blas_int* ldb,
          double*                        work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       info);
  void
  sormqr_(const char*                    side,
          const char*                    trans,
          const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          const dealii::types::blas_int* k,
          const float*                   A,
          const dealii::types::blas_int* lda,
          const float*                   tau,
          float*                         B,
          const dealii::types::blas_int* ldb,
          float*                         work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       info);
  // Compute matrix Q from the result of dgeqrf_
  void
  dorgqr_(const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          const dealii::types::blas_int* k,
          const double*                  A,
          const dealii::types::blas_int* lda,
          const double*                  tau,
          double*                        work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       info);
  void
  sorgqr_(const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          const dealii::types::blas_int* k,
          const float*                   A,
          const dealii::types::blas_int* lda,
          const float*                   tau,
          float*                         work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       info);
  // Compute Rx = b
  void
  dtrtrs_(const char*                    uplo,
          const char*                    trans,
          const char*                    diag,
          const dealii::types::blas_int* n,
          const dealii::types::blas_int* n_rhs,
          const double*                  A,
          const dealii::types::blas_int* lda,
          double*                        B,
          const dealii::types::blas_int* ldb,
          dealii::types::blas_int*       info);
  void
  strtrs_(const char*                    uplo,
          const char*                    trans,
          const char*                    diag,
          const dealii::types::blas_int* n,
          const dealii::types::blas_int* n_rhs,
          const float*                   A,
          const dealii::types::blas_int* lda,
          float*                         B,
          const dealii::types::blas_int* ldb,
          dealii::types::blas_int*       info);
  // Compute eigenvalues and vectors
  void
  dgeev_(const char*                    jobvl,
         const char*                    jobvr,
         const dealii::types::blas_int* n,
         double*                        A,
         const dealii::types::blas_int* lda,
         double*                        lambda_re,
         double*                        lambda_im,
         double*                        vl,
         const dealii::types::blas_int* ldvl,
         double*                        vr,
         const dealii::types::blas_int* ldva,
         double*                        work,
         const dealii::types::blas_int* lwork,
         dealii::types::blas_int*       info);
  void
  sgeev_(const char*                    jobvl,
         const char*                    jobvr,
         const dealii::types::blas_int* n,
         float*                         A,
         const dealii::types::blas_int* lda,
         float*                         lambda_re,
         float*                         lambda_im,
         float*                         vl,
         const dealii::types::blas_int* ldvl,
         float*                         vr,
         const dealii::types::blas_int* ldva,
         float*                         work,
         const dealii::types::blas_int* lwork,
         dealii::types::blas_int*       info);
  // Compute eigenvalues and vectors (expert)
  void
  dgeevx_(const char*                    balanc,
          const char*                    jobvl,
          const char*                    jobvr,
          const char*                    sense,
          const dealii::types::blas_int* n,
          double*                        A,
          const dealii::types::blas_int* lda,
          double*                        lambda_re,
          double*                        lambda_im,
          double*                        vl,
          const dealii::types::blas_int* ldvl,
          double*                        vr,
          const dealii::types::blas_int* ldvr,
          dealii::types::blas_int*       ilo,
          dealii::types::blas_int*       ihi,
          double*                        scale,
          double*                        abnrm,
          double*                        rconde,
          double*                        rcondv,
          double*                        work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       info);
  void
  sgeevx_(const char*                    balanc,
          const char*                    jobvl,
          const char*                    jobvr,
          const char*                    sense,
          const dealii::types::blas_int* n,
          float*                         A,
          const dealii::types::blas_int* lda,
          float*                         lambda_re,
          float*                         lambda_im,
          float*                         vl,
          const dealii::types::blas_int* ldvl,
          float*                         vr,
          const dealii::types::blas_int* ldvr,
          dealii::types::blas_int*       ilo,
          dealii::types::blas_int*       ihi,
          float*                         scale,
          float*                         abnrm,
          float*                         rconde,
          float*                         rcondv,
          float*                         work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       info);
  // Eigenvalues for a symmetric matrix
  void
  dsyev_(const char*                    jobz,
         const char*                    uplo,
         const dealii::types::blas_int* n,
         double*                        A,
         const dealii::types::blas_int* lda,
         double*                        w,
         double*                        work,
         const dealii::types::blas_int* lwork,
         dealii::types::blas_int*       info);
  void
  ssyev_(const char*                    jobz,
         const char*                    uplo,
         const dealii::types::blas_int* n,
         float*                         A,
         const dealii::types::blas_int* lda,
         float*                         w,
         float*                         work,
         const dealii::types::blas_int* lwork,
         dealii::types::blas_int*       info);
  // Same functionality as dsyev_ but with more options: E.g.
  // Compute only eigenvalues in a specific dealii::types::blas_interval,
  // Compute only eigenvalues with a specific index,
  // Set tolerance for eigenvalue computation
  void
  dsyevx_(const char*                    jobz,
          const char*                    range,
          const char*                    uplo,
          const dealii::types::blas_int* n,
          double*                        A,
          const dealii::types::blas_int* lda,
          const double*                  vl,
          const double*                  vu,
          const dealii::types::blas_int* il,
          const dealii::types::blas_int* iu,
          const double*                  abstol,
          dealii::types::blas_int*       m,
          double*                        w,
          double*                        z,
          const dealii::types::blas_int* ldz,
          double*                        work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       ifail,
          dealii::types::blas_int*       info);
  void
  ssyevx_(const char*                    jobz,
          const char*                    range,
          const char*                    uplo,
          const dealii::types::blas_int* n,
          float*                         A,
          const dealii::types::blas_int* lda,
          const float*                   vl,
          const float*                   vu,
          const dealii::types::blas_int* il,
          const dealii::types::blas_int* iu,
          const float*                   abstol,
          dealii::types::blas_int*       m,
          float*                         w,
          float*                         z,
          const dealii::types::blas_int* ldz,
          float*                         work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       ifail,
          dealii::types::blas_int*       info);
  // Same functionality as dsyev_ but using MRRR algorithm and with more options:
  // E.g. compute only eigenvalues in a specific dealii::types::blas_interval,
  // Compute only eigenvalues with a specific index.
  void
  dsyevr_(const char*                    jobz,
          const char*                    range,
          const char*                    uplo,
          const dealii::types::blas_int* n,
          double*                        A,
          const dealii::types::blas_int* lda,
          const double*                  vl,
          const double*                  vu,
          const dealii::types::blas_int* il,
          const dealii::types::blas_int* iu,
          const double*                  abstol,
          dealii::types::blas_int*       m,
          double*                        w,
          double*                        z,
          const dealii::types::blas_int* ldz,
          dealii::types::blas_int*       isuppz,
          double*                        work,
          dealii::types::blas_int*       lwork,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       liwork,
          dealii::types::blas_int*       info);
  void
  ssyevr_(const char*                    jobz,
          const char*                    range,
          const char*                    uplo,
          const dealii::types::blas_int* n,
          float*                         A,
          const dealii::types::blas_int* lda,
          const float*                   vl,
          const float*                   vu,
          const dealii::types::blas_int* il,
          const dealii::types::blas_int* iu,
          const float*                   abstol,
          dealii::types::blas_int*       m,
          float*                         w,
          float*                         z,
          const dealii::types::blas_int* ldz,
          dealii::types::blas_int*       isuppz,
          float*                         work,
          dealii::types::blas_int*       lwork,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       liwork,
          dealii::types::blas_int*       info);
  // Generalized eigenvalues and eigenvectors of
  // 1: A*x = lambda*B*x; 2: A*B*x = lambda*x; 3: B*A*x = lambda*x
  // A and B are symmetric and B is definite
  void
  dsygv_(const dealii::types::blas_int* itype,
         const char*                    jobz,
         const char*                    uplo,
         const dealii::types::blas_int* n,
         double*                        A,
         const dealii::types::blas_int* lda,
         double*                        B,
         const dealii::types::blas_int* ldb,
         double*                        w,
         double*                        work,
         const dealii::types::blas_int* lwork,
         dealii::types::blas_int*       info);
  void
  ssygv_(const dealii::types::blas_int* itype,
         const char*                    jobz,
         const char*                    uplo,
         const dealii::types::blas_int* n,
         float*                         A,
         const dealii::types::blas_int* lda,
         float*                         B,
         const dealii::types::blas_int* ldb,
         float*                         w,
         float*                         work,
         const dealii::types::blas_int* lwork,
         dealii::types::blas_int*       info);
  // Same functionality as dsygv_ but with more options: E.g.
  // Compute only eigenvalues in a specific dealii::types::blas_interval,
  // Compute only eigenvalues with a specific index,
  // Set tolerance for eigenvalue computation
  void
  dsygvx_(const dealii::types::blas_int* itype,
          const char*                    jobz,
          const char*                    range,
          const char*                    uplo,
          const dealii::types::blas_int* n,
          double*                        A,
          const dealii::types::blas_int* lda,
          double*                        B,
          const dealii::types::blas_int* ldb,
          const double*                  vl,
          const double*                  vu,
          const dealii::types::blas_int* il,
          const dealii::types::blas_int* iu,
          const double*                  abstol,
          dealii::types::blas_int*       m,
          double*                        w,
          double*                        z,
          const dealii::types::blas_int* ldz,
          double*                        work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       ifail,
          dealii::types::blas_int*       info);
  void
  ssygvx_(const dealii::types::blas_int* itype,
          const char*                    jobz,
          const char*                    range,
          const char*                    uplo,
          const dealii::types::blas_int* n,
          float*                         A,
          const dealii::types::blas_int* lda,
          float*                         B,
          const dealii::types::blas_int* ldb,
          const float*                   vl,
          const float*                   vu,
          const dealii::types::blas_int* il,
          const dealii::types::blas_int* iu,
          const float*                   abstol,
          dealii::types::blas_int*       m,
          float*                         w,
          float*                         z,
          const dealii::types::blas_int* ldz,
          float*                         work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       ifail,
          dealii::types::blas_int*       info);
  // Compute singular value decomposition using divide and conquer
  void
  dgesdd_(const char*                    jobz,
          const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          double*                        A,
          const dealii::types::blas_int* lda,
          double*                        s,
          double*                        u,
          const dealii::types::blas_int* ldu,
          double*                        vt,
          const dealii::types::blas_int* ldvt,
          double*                        work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       info);
  void
  sgesdd_(const char*                    jobz,
          const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          float*                         A,
          const dealii::types::blas_int* lda,
          float*                         s,
          float*                         u,
          const dealii::types::blas_int* ldu,
          float*                         vt,
          const dealii::types::blas_int* ldvt,
          float*                         work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       info);
  // Compute singular value decomposition
  void
  dgesvd_(dealii::types::blas_int*       jobu,
          dealii::types::blas_int*       jobvt,
          const dealii::types::blas_int* n,
          const dealii::types::blas_int* m,
          double*                        A,
          const dealii::types::blas_int* lda,
          double*                        s,
          double*                        u,
          const dealii::types::blas_int* ldu,
          double*                        vt,
          const dealii::types::blas_int* ldvt,
          double*                        work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       info);
  void
  sgesvd_(dealii::types::blas_int*       jobu,
          dealii::types::blas_int*       jobvt,
          const dealii::types::blas_int* n,
          const dealii::types::blas_int* m,
          float*                         A,
          const dealii::types::blas_int* lda,
          float*                         s,
          float*                         u,
          const dealii::types::blas_int* ldu,
          float*                         vt,
          const dealii::types::blas_int* ldvt,
          float*                         work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       info);
  // Solve a least squares problem using SVD
  void
  dgelsd_(const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          const dealii::types::blas_int* nrhs,
          const double*                  A,
          const dealii::types::blas_int* lda,
          double*                        B,
          const dealii::types::blas_int* ldb,
          double*                        s,
          const double*                  rcond,
          dealii::types::blas_int*       rank,
          double*                        work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       info);
  void
  sgelsd_(const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          const dealii::types::blas_int* nrhs,
          const float*                   A,
          const dealii::types::blas_int* lda,
          float*                         B,
          const dealii::types::blas_int* ldb,
          float*                         s,
          const float*                   rcond,
          dealii::types::blas_int*       rank,
          float*                         work,
          const dealii::types::blas_int* lwork,
          dealii::types::blas_int*       iwork,
          dealii::types::blas_int*       info);
  // Symmetric tridiagonal matrix
  void
  dstev_(const char*                    jobz,
         const dealii::types::blas_int* n,
         double*                        d,
         double*                        e,
         double*                        z,
         const dealii::types::blas_int* ldz,
         double*                        work,
         dealii::types::blas_int*       info);
  void
  sstev_(const char*                    jobz,
         const dealii::types::blas_int* n,
         float*                         d,
         float*                         e,
         float*                         z,
         const dealii::types::blas_int* ldz,
         float*                         work,
         dealii::types::blas_int*       info);
  // Rank-1 update for symmetric matrices
  void
  dsyr_(const char*                    uplo,
        const dealii::types::blas_int* n,
        const double*                  alpha,
        const double*                  x,
        const dealii::types::blas_int* incx,
        double*                        A,
        const dealii::types::blas_int* lda);
  void
  ssyr_(const char*                    uplo,
        const dealii::types::blas_int* n,
        const float*                   alpha,
        const float*                   x,
        const dealii::types::blas_int* incx,
        float*                         A,
        const dealii::types::blas_int* lda);
  // Multiply rectangular mxn real matrix by real scalar CTO/CFROM
  void
  dlascl_(const char*                    type,
          const dealii::types::blas_int* kl,
          const dealii::types::blas_int* ku,
          const double*                  cfrom,
          const double*                  cto,
          const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          double*                        A,
          const dealii::types::blas_int* lda,
          dealii::types::blas_int*       info);
  void
  slascl_(const char*                    type,
          const dealii::types::blas_int* kl,
          const dealii::types::blas_int* ku,
          const float*                   cfrom,
          const float*                   cto,
          const dealii::types::blas_int* m,
          const dealii::types::blas_int* n,
          float*                         A,
          const dealii::types::blas_int* lda,
          dealii::types::blas_int*       info);
  // dlamch and slamch help determining machine precision
  double
  dlamch_(const char* chmach);
  float
  slamch_(const char* chmach);
}

DEAL_II_NAMESPACE_OPEN

/// Template wrapper for LAPACK functions dsyrk and ssyrk
template <typename number>
inline void
syrk(const char*,
     const char*,
     const types::blas_int*,
     const types::blas_int*,
     const number*,
     const number*,
     const types::blas_int*,
     const number*,
     number*,
     const types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
syrk(const char*            uplo,
     const char*            trans,
     const types::blas_int* n,
     const types::blas_int* k,
     const double*          alpha,
     const double*          A,
     const types::blas_int* lda,
     const double*          beta,
     double*                C,
     const types::blas_int* ldc)
{
  dsyrk_(uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}

inline void
syrk(const char*            uplo,
     const char*            trans,
     const types::blas_int* n,
     const types::blas_int* k,
     const float*           alpha,
     const float*           A,
     const types::blas_int* lda,
     const float*           beta,
     float*                 C,
     const types::blas_int* ldc)
{
  ssyrk_(uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
#endif

/// Template wrapper for LAPACK functions dsyr and ssyr
template <typename number>
inline void
syr(const char*,
    const types::blas_int*,
    const number*,
    const number*,
    const types::blas_int*,
    number*,
    const types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
syr(const char*            uplo,
    const types::blas_int* n,
    const double*          alpha,
    const double*          x,
    const types::blas_int* incx,
    double*                A,
    const types::blas_int* lda)
{
  dsyr_(uplo, n, alpha, x, incx, A, lda);
}
inline void
syr(const char*            uplo,
    const types::blas_int* n,
    const float*           alpha,
    const float*           x,
    const types::blas_int* incx,
    float*                 A,
    const types::blas_int* lda)
{
  ssyr_(uplo, n, alpha, x, incx, A, lda);
}
#endif

/// Template wrapper for LAPACK functions daxpy and saxpy
template <typename number1, typename number2, typename number3>
inline void
axpy(const types::blas_int*,
     const number1*,
     const number2*,
     const types::blas_int*,
     number3*,
     const types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
axpy(const types::blas_int* n,
     const double*          alpha,
     const double*          x,
     const types::blas_int* incx,
     double*                y,
     const types::blas_int* incy)
{
  daxpy_(n, alpha, x, incx, y, incy);
}
#else
inline void
axpy(const types::blas_int*,
     const double*,
     const double*,
     const types::blas_int*,
     double*,
     const types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("daxpy"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
axpy(const types::blas_int* n,
     const float*           alpha,
     const float*           x,
     const types::blas_int* incx,
     float*                 y,
     const types::blas_int* incy)
{
  saxpy_(n, alpha, x, incx, y, incy);
}
#else
inline void
axpy(const types::blas_int*,
     const float*,
     const float*,
     const types::blas_int*,
     float*,
     const types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("saxpy"));
}
#endif

/// Template wrapper for LAPACK functions dgemv and sgemv
template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5>
inline void
gemv(const char*,
     const types::blas_int*,
     const types::blas_int*,
     const number1*,
     const number2*,
     const types::blas_int*,
     const number3*,
     const types::blas_int*,
     const number4*,
     number5*,
     const types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
gemv(const char*            trans,
     const types::blas_int* m,
     const types::blas_int* n,
     const double*          alpha,
     const double*          A,
     const types::blas_int* lda,
     const double*          x,
     const types::blas_int* incx,
     const double*          b,
     double*                y,
     const types::blas_int* incy)
{
  dgemv_(trans, m, n, alpha, A, lda, x, incx, b, y, incy);
}
#else
inline void
gemv(const char*,
     const types::blas_int*,
     const types::blas_int*,
     const double*,
     const double*,
     const types::blas_int*,
     const double*,
     const types::blas_int*,
     const double*,
     double*,
     const types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dgemv"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
gemv(const char*            trans,
     const types::blas_int* m,
     const types::blas_int* n,
     const float*           alpha,
     const float*           A,
     const types::blas_int* lda,
     const float*           x,
     const types::blas_int* incx,
     const float*           b,
     float*                 y,
     const types::blas_int* incy)
{
  sgemv_(trans, m, n, alpha, A, lda, x, incx, b, y, incy);
}
#else
inline void
gemv(const char*,
     const types::blas_int*,
     const types::blas_int*,
     const float*,
     const float*,
     const types::blas_int*,
     const float*,
     const types::blas_int*,
     const float*,
     float*,
     const types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("sgemv"));
}
#endif

/// Template wrapper for LAPACK functions dtrmv and strmv
template <typename number>
inline void
trmv(const char* /*uplo*/,
     const char* /*trans*/,
     const char* /*diag*/,
     const types::blas_int* /*N*/,
     const number* /*A*/,
     const types::blas_int* /*lda*/,
     number* /*x*/,
     const types::blas_int* /*incx*/)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
trmv(const char*            uplo,
     const char*            trans,
     const char*            diag,
     const types::blas_int* N,
     const double*          A,
     const types::blas_int* lda,
     double*                x,
     const types::blas_int* incx)
{
  dtrmv_(uplo, trans, diag, N, A, lda, x, incx);
}
#else
inline void
trmv(const char* /*uplo*/,
     const char* /*trans*/,
     const char* /*diag*/,
     const types::blas_int* /*N*/,
     const double* /*A*/,
     const types::blas_int* /*lda*/,
     double* /*x*/,
     const types::blas_int* /*incx*/)
{
  Assert(false, LAPACKSupport::ExcMissing("dtrmv"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
trmv(const char*            uplo,
     const char*            trans,
     const char*            diag,
     const types::blas_int* N,
     const float*           A,
     const types::blas_int* lda,
     float*                 x,
     const types::blas_int* incx)
{
  strmv_(uplo, trans, diag, N, A, lda, x, incx);
}
#else
inline void
trmv(const char* /*uplo*/,
     const char* /*trans*/,
     const char* /*diag*/,
     const types::blas_int* /*N*/,
     const float* /*A*/,
     const types::blas_int* /*lda*/,
     float* /*x*/,
     const types::blas_int* /*incx*/)
{
  Assert(false, LAPACKSupport::ExcMissing("dtrmv"));
}
#endif

/// Template wrapper for LAPACK functions dgemm and sgemm
template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5>
inline void
gemm(const char*,
     const char*,
     const types::blas_int*,
     const types::blas_int*,
     const types::blas_int*,
     const number1*,
     const number2*,
     const types::blas_int*,
     const number3*,
     const types::blas_int*,
     const number4*,
     number5*,
     const types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
gemm(const char*            transa,
     const char*            transb,
     const types::blas_int* m,
     const types::blas_int* n,
     const types::blas_int* k,
     const double*          alpha,
     const double*          A,
     const types::blas_int* lda,
     const double*          B,
     const types::blas_int* ldb,
     const double*          beta,
     double*                C,
     const types::blas_int* ldc)
{
  dgemm_(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
#else
inline void
gemm(const char*,
     const char*,
     const types::blas_int*,
     const types::blas_int*,
     const types::blas_int*,
     const double*,
     const double*,
     const types::blas_int*,
     const double*,
     const types::blas_int*,
     const double*,
     double*,
     const types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dgemm"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
gemm(const char*            transa,
     const char*            transb,
     const types::blas_int* m,
     const types::blas_int* n,
     const types::blas_int* k,
     const float*           alpha,
     const float*           A,
     const types::blas_int* lda,
     const float*           B,
     const types::blas_int* ldb,
     const float*           beta,
     float*                 C,
     const types::blas_int* ldc)
{
  sgemm_(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
#else
inline void
gemm(const char*,
     const char*,
     const types::blas_int*,
     const types::blas_int*,
     const types::blas_int*,
     const float*,
     const float*,
     const types::blas_int*,
     const float*,
     const types::blas_int*,
     const float*,
     float*,
     const types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("sgemm"));
}
#endif

/// Template wrapper for potrf
template <typename number1>
inline void
potrf(const char*,
      const types::blas_int*,
      number1*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

inline void
potrf(const char*            uplo,
      const types::blas_int* n,
      double*                A,
      const types::blas_int* lda,
      types::blas_int*       info)
{
#ifdef DEAL_II_WITH_LAPACK
  dpotrf_(uplo, n, A, lda, info);
#else
  (void) uplo;
  (void) n;
  (void) A;
  (void) lda;
  (void) info;

  Assert(false, LAPACKSupport::ExcMissing("dpotrf"));
#endif
}

inline void
potrf(const char*            uplo,
      const types::blas_int* n,
      float*                 A,
      const types::blas_int* lda,
      types::blas_int*       info)
{
#ifdef DEAL_II_WITH_LAPACK
  spotrf_(uplo, n, A, lda, info);
#else
  (void) uplo;
  (void) n;
  (void) A;
  (void) lda;
  (void) info;

  Assert(false, LAPACKSupport::ExcMissing("spotrf"));
#endif
}

/// Template wrapper for trcon
template <typename number>
inline void
trcon(const char* /*norm*/,
      const char* /*uplo*/,
      const char* /*diag*/,
      const types::blas_int* /*n*/,
      const number* /*A*/,
      const types::blas_int* /*lda*/,
      number* /*rcond*/,
      number* /*work*/,
      types::blas_int* /*iwork*/,
      types::blas_int* /*info*/)
{
  Assert(false, ExcNotImplemented());
}

inline void
trcon(const char*            norm,
      const char*            uplo,
      const char*            diag,
      const types::blas_int* n,
      const double*          A,
      const types::blas_int* lda,
      double*                rcond,
      double*                work,
      types::blas_int*       iwork,
      types::blas_int*       info)
{
#ifdef DEAL_II_WITH_LAPACK
  dtrcon_(norm, uplo, diag, n, A, lda, rcond, work, iwork, info);
#else
  (void) norm;
  (void) uplo;
  (void) diag;
  (void) n;
  (void) A;
  (void) lda;
  (void) rcond;
  (void) work;
  (void) iwork;
  (void) info;

  Assert(false, LAPACKSupport::ExcMissing("dtrcon"));
#endif
}

inline void
trcon(const char*            norm,
      const char*            uplo,
      const char*            diag,
      const types::blas_int* n,
      const float*           A,
      const types::blas_int* lda,
      float*                 rcond,
      float*                 work,
      types::blas_int*       iwork,
      types::blas_int*       info)
{
#ifdef DEAL_II_WITH_LAPACK
  strcon_(norm, uplo, diag, n, A, lda, rcond, work, iwork, info);
#else
  (void) norm;
  (void) uplo;
  (void) diag;
  (void) n;
  (void) A;
  (void) lda;
  (void) rcond;
  (void) work;
  (void) iwork;
  (void) info;

  Assert(false, LAPACKSupport::ExcMissing("dtrcon"));
#endif
}

/// Template wrapper for pocon
template <typename number1>
inline void
pocon(const char*,
      const types::blas_int*,
      const number1*,
      const types::blas_int*,
      const number1*,
      number1*,
      number1*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

inline void
pocon(const char*            uplo,
      const types::blas_int* n,
      const double*          A,
      const types::blas_int* lda,
      const double*          anorm,
      double*                rcond,
      double*                work,
      types::blas_int*       iwork,
      types::blas_int*       info)
{
#ifdef DEAL_II_WITH_LAPACK
  dpocon_(uplo, n, A, lda, anorm, rcond, work, iwork, info);
#else
  (void) uplo;
  (void) n;
  (void) A;
  (void) lda;
  (void) anorm;
  (void) rcond;
  (void) work;
  (void) iwork;
  (void) info;

  Assert(false, LAPACKSupport::ExcMissing("dpocon"));
#endif
}

inline void
pocon(const char*            uplo,
      const types::blas_int* n,
      const float*           A,
      const types::blas_int* lda,
      const float*           anorm,
      float*                 rcond,
      float*                 work,
      types::blas_int*       iwork,
      types::blas_int*       info)
{
#ifdef DEAL_II_WITH_LAPACK
  spocon_(uplo, n, A, lda, anorm, rcond, work, iwork, info);
#else
  (void) uplo;
  (void) n;
  (void) A;
  (void) lda;
  (void) anorm;
  (void) rcond;
  (void) work;
  (void) iwork;
  (void) info;

  Assert(false, LAPACKSupport::ExcMissing("dpocon"));
#endif
}

/// Template wrapper for potri
template <typename number1>
inline void
potri(const char*,
      const types::blas_int*,
      number1*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

inline void
potri(const char*            uplo,
      const types::blas_int* n,
      double*                A,
      const types::blas_int* lda,
      types::blas_int*       info)
{
#ifdef DEAL_II_WITH_LAPACK
  dpotri_(uplo, n, A, lda, info);
#else
  (void) uplo;
  (void) n;
  (void) A;
  (void) lda;
  (void) info;

  Assert(false, LAPACKSupport::ExcMissing("dpotri"));
#endif
}

inline void
potri(const char*            uplo,
      const types::blas_int* n,
      float*                 A,
      const types::blas_int* lda,
      types::blas_int*       info)
{
#ifdef DEAL_II_WITH_LAPACK
  spotri_(uplo, n, A, lda, info);
#else
  (void) uplo;
  (void) n;
  (void) A;
  (void) lda;
  (void) info;

  Assert(false, LAPACKSupport::ExcMissing("spotri"));
#endif
}

/// Template wrapper for lansy
template <typename number>
inline number
lansy(const char*,
      const char*,
      const types::blas_int*,
      const number*,
      const types::blas_int*,
      number*)
{
  Assert(false, ExcNotImplemented());
  return number();
}

inline double
lansy(const char*            norm,
      const char*            uplo,
      const types::blas_int* n,
      const double*          A,
      const types::blas_int* lda,
      double*                work)
{
#ifdef DEAL_II_WITH_LAPACK
  return dlansy_(norm, uplo, n, A, lda, work);
#else
  (void) norm;
  (void) uplo;
  (void) n;
  (void) A;
  (void) lda;
  (void) work;

  Assert(false, LAPACKSupport::ExcMissing("lansy"));
  return 0.;
#endif
}

inline float
lansy(const char*            norm,
      const char*            uplo,
      const types::blas_int* n,
      const float*           A,
      const types::blas_int* lda,
      float*                 work)
{
#ifdef DEAL_II_WITH_LAPACK
  return slansy_(norm, uplo, n, A, lda, work);
#else
  (void) norm;
  (void) uplo;
  (void) n;
  (void) A;
  (void) lda;
  (void) work;

  Assert(false, LAPACKSupport::ExcMissing("lansy"));
  return 0.;
#endif
}

/// Template wrapper for lange
template <typename number>
inline number
lange(const char*,
      const types::blas_int*,
      const types::blas_int*,
      const number*,
      const types::blas_int*,
      number*)
{
  Assert(false, ExcNotImplemented());
  return number();
}

inline double
lange(const char*            norm,
      const types::blas_int* m,
      const types::blas_int* n,
      const double*          A,
      const types::blas_int* lda,
      double*                work)
{
#ifdef DEAL_II_WITH_LAPACK
  return dlange_(norm, m, n, A, lda, work);
#else
  (void) norm;
  (void) m;
  (void) n;
  (void) A;
  (void) lda;
  (void) work;

  Assert(false, LAPACKSupport::ExcMissing("lange"));
  return 0.;
#endif
}

inline float
lange(const char*            norm,
      const types::blas_int* m,
      const types::blas_int* n,
      const float*           A,
      const types::blas_int* lda,
      float*                 work)
{
#ifdef DEAL_II_WITH_LAPACK
  return slange_(norm, m, n, A, lda, work);
#else
  (void) norm;
  (void) m;
  (void) n;
  (void) A;
  (void) lda;
  (void) work;

  Assert(false, LAPACKSupport::ExcMissing("lange"));
  return 0.;
#endif
}

/// Template wrapper for LAPACK functions dgetrf and sgetrf
template <typename number1>
inline void
getrf(const types::blas_int*,
      const types::blas_int*,
      number1*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
getrf(const types::blas_int* m,
      const types::blas_int* n,
      double*                A,
      const types::blas_int* lda,
      types::blas_int*       ipiv,
      types::blas_int*       info)
{
  dgetrf_(m, n, A, lda, ipiv, info);
}
#else
inline void
getrf(const types::blas_int*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dgetrf"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
getrf(const types::blas_int* m,
      const types::blas_int* n,
      float*                 A,
      const types::blas_int* lda,
      types::blas_int*       ipiv,
      types::blas_int*       info)
{
  sgetrf_(m, n, A, lda, ipiv, info);
}
#else
inline void
getrf(const types::blas_int*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("sgetrf"));
}
#endif

/// Template wrapper for LAPACK functions dgetrs and sgetrs
template <typename number1, typename number2>
inline void
getrs(const char*,
      const types::blas_int*,
      const types::blas_int*,
      const number1*,
      const types::blas_int*,
      const types::blas_int*,
      number2*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
getrs(const char*            trans,
      const types::blas_int* n,
      const types::blas_int* nrhs,
      const double*          A,
      const types::blas_int* lda,
      const types::blas_int* ipiv,
      double*                b,
      const types::blas_int* ldb,
      types::blas_int*       info)
{
  dgetrs_(trans, n, nrhs, A, lda, ipiv, b, ldb, info);
}
#else
inline void
getrs(const char*,
      const types::blas_int*,
      const types::blas_int*,
      const double*,
      const types::blas_int*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dgetrs"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
getrs(const char*            trans,
      const types::blas_int* n,
      const types::blas_int* nrhs,
      const float*           A,
      const types::blas_int* lda,
      const types::blas_int* ipiv,
      float*                 b,
      const types::blas_int* ldb,
      types::blas_int*       info)
{
  sgetrs_(trans, n, nrhs, A, lda, ipiv, b, ldb, info);
}
#else
inline void
getrs(const char*,
      const types::blas_int*,
      const types::blas_int*,
      const float*,
      const types::blas_int*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("sgetrs"));
}
#endif

///  Template wrapper for LAPACK functions dpotrs and spotrs
template <typename number>
inline void
potrs(const char*,
      const types::blas_int*,
      const types::blas_int*,
      const number*,
      const types::blas_int*,
      number*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
potrs(const char*            uplo,
      const types::blas_int* n,
      const types::blas_int* nrhs,
      const double*          A,
      const types::blas_int* lda,
      double*                B,
      const types::blas_int* ldb,
      types::blas_int*       info)
{
  dpotrs_(uplo, n, nrhs, A, lda, B, ldb, info);
}
inline void
potrs(const char*            uplo,
      const types::blas_int* n,
      const types::blas_int* nrhs,
      const float*           A,
      const types::blas_int* lda,
      float*                 B,
      const types::blas_int* ldb,
      types::blas_int*       info)
{
  spotrs_(uplo, n, nrhs, A, lda, B, ldb, info);
}
#else
inline void
potrs(const char*,
      const types::blas_int*,
      const types::blas_int*,
      const double*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dpotrs"));
}
inline void
potrs(const char*,
      const types::blas_int*,
      const types::blas_int*,
      const float*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("spotrs"));
}
#endif

/// Template wrapper for LAPACK functions dgetri and sgetri
template <typename number1, typename number2>
inline void
getri(const types::blas_int*,
      number1*,
      const types::blas_int*,
      types::blas_int*,
      number2*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
getri(const types::blas_int* n,
      double*                A,
      const types::blas_int* lda,
      types::blas_int*       ipiv,
      double*                inv_work,
      const types::blas_int* lwork,
      types::blas_int*       info)
{
  dgetri_(n, A, lda, ipiv, inv_work, lwork, info);
}
#else
inline void
getri(const types::blas_int*,
      double*,
      const types::blas_int*,
      types::blas_int*,
      double*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dgetri"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
getri(const types::blas_int* n,
      float*                 A,
      const types::blas_int* lda,
      types::blas_int*       ipiv,
      float*                 inv_work,
      const types::blas_int* lwork,
      types::blas_int*       info)
{
  sgetri_(n, A, lda, ipiv, inv_work, lwork, info);
}
#else
inline void
getri(const types::blas_int*,
      float*,
      const types::blas_int*,
      types::blas_int*,
      float*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("sgetri"));
}
#endif

/// Template wrapper for LAPACK functions dgeqrf and sgeqrf
template <typename number1, typename number2, typename number3>
inline void
geqrf(const types::blas_int*,
      const types::blas_int*,
      number1*,
      const types::blas_int*,
      number2*,
      number3*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
geqrf(const types::blas_int* m,
      const types::blas_int* n,
      double*                A,
      const types::blas_int* lda,
      double*                tau,
      double*                work,
      const types::blas_int* lwork,
      types::blas_int*       info)
{
  dgeqrf_(m, n, A, lda, tau, work, lwork, info);
}
#else
inline void
geqrf(const types::blas_int*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      double*,
      double*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dgeqrf"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
geqrf(const types::blas_int* m,
      const types::blas_int* n,
      float*                 A,
      const types::blas_int* lda,
      float*                 tau,
      float*                 work,
      const types::blas_int* lwork,
      types::blas_int*       info)
{
  sgeqrf_(m, n, A, lda, tau, work, lwork, info);
}
#else
inline void
geqrf(const types::blas_int*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      float*,
      float*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("sgeqrf"));
}
#endif

/// Template wrapper for LAPACK functions dormqr and sormqr
template <typename number1,
          typename number2,
          typename number3,
          typename number4>
inline void
ormqr(const char*,
      const char*,
      const types::blas_int*,
      const types::blas_int*,
      const types::blas_int*,
      const number1*,
      const types::blas_int*,
      const number2*,
      number3*,
      const types::blas_int*,
      number4*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
ormqr(const char*            side,
      const char*            trans,
      const types::blas_int* m,
      const types::blas_int* n,
      const types::blas_int* k,
      const double*          A,
      const types::blas_int* lda,
      const double*          tau,
      double*                B,
      const types::blas_int* ldb,
      double*                work,
      const types::blas_int* lwork,
      types::blas_int*       info)
{
  dormqr_(side, trans, m, n, k, A, lda, tau, B, ldb, work, lwork, info);
}
#else
inline void
ormqr(const char*,
      const char*,
      const types::blas_int*,
      const types::blas_int*,
      const types::blas_int*,
      const double*,
      const types::blas_int*,
      const double*,
      double*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dormqr"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
ormqr(const char*            side,
      const char*            trans,
      const types::blas_int* m,
      const types::blas_int* n,
      const types::blas_int* k,
      const float*           A,
      const types::blas_int* lda,
      const float*           tau,
      float*                 B,
      const types::blas_int* ldb,
      float*                 work,
      const types::blas_int* lwork,
      types::blas_int*       info)
{
  sormqr_(side, trans, m, n, k, A, lda, tau, B, ldb, work, lwork, info);
}
#else
inline void
ormqr(const char*,
      const char*,
      const types::blas_int*,
      const types::blas_int*,
      const types::blas_int*,
      const float*,
      const types::blas_int*,
      const float*,
      float*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("sormqr"));
}
#endif

/// Template wrapper for LAPACK functions dorgqr and sorgqr
template <typename number1, typename number2, typename number3>
inline void
orgqr(const types::blas_int*,
      const types::blas_int*,
      const types::blas_int*,
      const number1*,
      const types::blas_int*,
      const number2*,
      number3*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
orgqr(const types::blas_int* m,
      const types::blas_int* n,
      const types::blas_int* k,
      const double*          A,
      const types::blas_int* lda,
      const double*          tau,
      double*                work,
      const types::blas_int* lwork,
      types::blas_int*       info)
{
  dorgqr_(m, n, k, A, lda, tau, work, lwork, info);
}
#else
inline void
orgqr(const types::blas_int*,
      const types::blas_int*,
      const types::blas_int*,
      const double*,
      const types::blas_int*,
      const double*,
      double*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dorgqr"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
orgqr(const types::blas_int* m,
      const types::blas_int* n,
      const types::blas_int* k,
      const float*           A,
      const types::blas_int* lda,
      const float*           tau,
      float*                 work,
      const types::blas_int* lwork,
      types::blas_int*       info)
{
  sorgqr_(m, n, k, A, lda, tau, work, lwork, info);
}
#else
inline void
orgqr(const types::blas_int*,
      const types::blas_int*,
      const types::blas_int*,
      const float*,
      const types::blas_int*,
      const float*,
      float*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("sorgqr"));
}
#endif

/// Template wrapper for LAPACK functions dtrtrs and strtrs
template <typename number1, typename number2>
inline void
trtrs(const char*,
      const char*,
      const char*,
      const types::blas_int*,
      const types::blas_int*,
      const number1*,
      const types::blas_int*,
      number2*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
trtrs(const char*            uplo,
      const char*            trans,
      const char*            diag,
      const types::blas_int* n,
      const types::blas_int* n_rhs,
      const double*          A,
      const types::blas_int* lda,
      double*                B,
      const types::blas_int* ldb,
      types::blas_int*       info)
{
  dtrtrs_(uplo, trans, diag, n, n_rhs, A, lda, B, ldb, info);
}
#else
inline void
trtrs(const char*,
      const char*,
      const char*,
      const types::blas_int*,
      const types::blas_int*,
      const double*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dtrtrs"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
trtrs(const char*            uplo,
      const char*            trans,
      const char*            diag,
      const types::blas_int* n,
      const types::blas_int* n_rhs,
      const float*           A,
      const types::blas_int* lda,
      float*                 B,
      const types::blas_int* ldb,
      types::blas_int*       info)
{
  strtrs_(uplo, trans, diag, n, n_rhs, A, lda, B, ldb, info);
}
#else
inline void
trtrs(const char*,
      const char*,
      const char*,
      const types::blas_int*,
      const types::blas_int*,
      const float*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("strtrs"));
}
#endif

/// Template wrapper for LAPACK functions dgeev and sgeev
template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5,
          typename number6>
inline void
geev(const char*,
     const char*,
     const types::blas_int*,
     number1*,
     const types::blas_int*,
     number2*,
     number3*,
     number4*,
     const types::blas_int*,
     number5*,
     const types::blas_int*,
     number6*,
     const types::blas_int*,
     types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
geev(const char*            jobvl,
     const char*            jobvr,
     const types::blas_int* n,
     double*                A,
     const types::blas_int* lda,
     double*                lambda_re,
     double*                lambda_im,
     double*                vl,
     const types::blas_int* ldvl,
     double*                vr,
     const types::blas_int* ldva,
     double*                work,
     const types::blas_int* lwork,
     types::blas_int*       info)
{
  dgeev_(jobvl,
         jobvr,
         n,
         A,
         lda,
         lambda_re,
         lambda_im,
         vl,
         ldvl,
         vr,
         ldva,
         work,
         lwork,
         info);
}
#else
inline void
geev(const char*,
     const char*,
     const types::blas_int*,
     double*,
     const types::blas_int*,
     double*,
     double*,
     double*,
     const types::blas_int*,
     double*,
     const types::blas_int*,
     double*,
     const types::blas_int*,
     types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dgeev"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
geev(const char*            jobvl,
     const char*            jobvr,
     const types::blas_int* n,
     float*                 A,
     const types::blas_int* lda,
     float*                 lambda_re,
     float*                 lambda_im,
     float*                 vl,
     const types::blas_int* ldvl,
     float*                 vr,
     const types::blas_int* ldva,
     float*                 work,
     const types::blas_int* lwork,
     types::blas_int*       info)
{
  sgeev_(jobvl,
         jobvr,
         n,
         A,
         lda,
         lambda_re,
         lambda_im,
         vl,
         ldvl,
         vr,
         ldva,
         work,
         lwork,
         info);
}
#else
inline void
geev(const char*,
     const char*,
     const types::blas_int*,
     float*,
     const types::blas_int*,
     float*,
     float*,
     float*,
     const types::blas_int*,
     float*,
     const types::blas_int*,
     float*,
     const types::blas_int*,
     types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("sgeev"));
}
#endif

/// Template wrapper for LAPACK functions dgeevx and sgeevx
template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5,
          typename number6,
          typename number7,
          typename number8,
          typename number9,
          typename number10>
inline void
geevx(const char*,
      const char*,
      const char*,
      const char*,
      const types::blas_int*,
      number1*,
      const types::blas_int*,
      number2*,
      number3*,
      number4*,
      const types::blas_int*,
      number5*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*,
      number6*,
      number7*,
      number8*,
      number9*,
      number10*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
geevx(const char*            balanc,
      const char*            jobvl,
      const char*            jobvr,
      const char*            sense,
      const types::blas_int* n,
      double*                A,
      const types::blas_int* lda,
      double*                lambda_re,
      double*                lambda_im,
      double*                vl,
      const types::blas_int* ldvl,
      double*                vr,
      const types::blas_int* ldvr,
      types::blas_int*       ilo,
      types::blas_int*       ihi,
      double*                scale,
      double*                abnrm,
      double*                rconde,
      double*                rcondv,
      double*                work,
      const types::blas_int* lwork,
      types::blas_int*       iwork,
      types::blas_int*       info)
{
  dgeevx_(balanc,
          jobvl,
          jobvr,
          sense,
          n,
          A,
          lda,
          lambda_re,
          lambda_im,
          vl,
          ldvl,
          vr,
          ldvr,
          ilo,
          ihi,
          scale,
          abnrm,
          rconde,
          rcondv,
          work,
          lwork,
          iwork,
          info);
}
#else
inline void
geevx(const char*,
      const char*,
      const char*,
      const char*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      double*,
      double*,
      double*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*,
      double*,
      double*,
      double*,
      double*,
      double*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dgeevx"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
geevx(const char*            balanc,
      const char*            jobvl,
      const char*            jobvr,
      const char*            sense,
      const types::blas_int* n,
      float*                 A,
      const types::blas_int* lda,
      float*                 lambda_re,
      float*                 lambda_im,
      float*                 vl,
      const types::blas_int* ldvl,
      float*                 vr,
      const types::blas_int* ldvr,
      types::blas_int*       ilo,
      types::blas_int*       ihi,
      float*                 scale,
      float*                 abnrm,
      float*                 rconde,
      float*                 rcondv,
      float*                 work,
      const types::blas_int* lwork,
      types::blas_int*       iwork,
      types::blas_int*       info)
{
  sgeevx_(balanc,
          jobvl,
          jobvr,
          sense,
          n,
          A,
          lda,
          lambda_re,
          lambda_im,
          vl,
          ldvl,
          vr,
          ldvr,
          ilo,
          ihi,
          scale,
          abnrm,
          rconde,
          rcondv,
          work,
          lwork,
          iwork,
          info);
}
#else
inline void
geevx(const char*,
      const char*,
      const char*,
      const char*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      float*,
      float*,
      float*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*,
      float*,
      float*,
      float*,
      float*,
      float*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("sgeevx"));
}
#endif

/// Template wrapper for LAPACK functions dsyev and ssyev
template <typename number1, typename number2, typename number3>
inline void
syev(const char*,
     const char*,
     const types::blas_int*,
     number1*,
     const types::blas_int*,
     number2*,
     number3*,
     const types::blas_int*,
     types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
syev(const char*            jobz,
     const char*            uplo,
     const types::blas_int* n,
     double*                A,
     const types::blas_int* lda,
     double*                w,
     double*                work,
     const types::blas_int* lwork,
     types::blas_int*       info)
{
  dsyev_(jobz, uplo, n, A, lda, w, work, lwork, info);
}
#else
inline void
syev(const char*,
     const char*,
     const types::blas_int*,
     double*,
     const types::blas_int*,
     double*,
     double*,
     const types::blas_int*,
     types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dsyev"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
syev(const char*            jobz,
     const char*            uplo,
     const types::blas_int* n,
     float*                 A,
     const types::blas_int* lda,
     float*                 w,
     float*                 work,
     const types::blas_int* lwork,
     types::blas_int*       info)
{
  ssyev_(jobz, uplo, n, A, lda, w, work, lwork, info);
}
#else
inline void
syev(const char*,
     const char*,
     const types::blas_int*,
     float*,
     const types::blas_int*,
     float*,
     float*,
     const types::blas_int*,
     types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("ssyev"));
}
#endif

/// Template wrapper for LAPACK functions dsyevx and ssyevx
template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5,
          typename number6,
          typename number7>
inline void
syevx(const char*,
      const char*,
      const char*,
      const types::blas_int*,
      number1*,
      const types::blas_int*,
      const number2*,
      const number3*,
      const types::blas_int*,
      const types::blas_int*,
      const number4*,
      types::blas_int*,
      number5*,
      number6*,
      const types::blas_int*,
      number7*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
syevx(const char*            jobz,
      const char*            range,
      const char*            uplo,
      const types::blas_int* n,
      double*                A,
      const types::blas_int* lda,
      const double*          vl,
      const double*          vu,
      const types::blas_int* il,
      const types::blas_int* iu,
      const double*          abstol,
      types::blas_int*       m,
      double*                w,
      double*                z,
      const types::blas_int* ldz,
      double*                work,
      const types::blas_int* lwork,
      types::blas_int*       iwork,
      types::blas_int*       ifail,
      types::blas_int*       info)
{
  dsyevx_(jobz,
          range,
          uplo,
          n,
          A,
          lda,
          vl,
          vu,
          il,
          iu,
          abstol,
          m,
          w,
          z,
          ldz,
          work,
          lwork,
          iwork,
          ifail,
          info);
}
#else
inline void
syevx(const char*,
      const char*,
      const char*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      const double*,
      const double*,
      const types::blas_int*,
      const types::blas_int*,
      const double*,
      types::blas_int*,
      double*,
      double*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dsyevx"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
syevx(const char*            jobz,
      const char*            range,
      const char*            uplo,
      const types::blas_int* n,
      float*                 A,
      const types::blas_int* lda,
      const float*           vl,
      const float*           vu,
      const types::blas_int* il,
      const types::blas_int* iu,
      const float*           abstol,
      types::blas_int*       m,
      float*                 w,
      float*                 z,
      const types::blas_int* ldz,
      float*                 work,
      const types::blas_int* lwork,
      types::blas_int*       iwork,
      types::blas_int*       ifail,
      types::blas_int*       info)
{
  ssyevx_(jobz,
          range,
          uplo,
          n,
          A,
          lda,
          vl,
          vu,
          il,
          iu,
          abstol,
          m,
          w,
          z,
          ldz,
          work,
          lwork,
          iwork,
          ifail,
          info);
}
#else
inline void
syevx(const char*,
      const char*,
      const char*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      const float*,
      const float*,
      const types::blas_int*,
      const types::blas_int*,
      const float*,
      types::blas_int*,
      float*,
      float*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("ssyevx"));
}
#endif

// Template wrapper for LAPACK functions dsyevr and ssyevr
template <typename number>
inline void
syevr(const char* /*jobz*/,
      const char* /*range*/,
      const char* /*uplo*/,
      const types::blas_int* /*n*/,
      number* /*A*/,
      const types::blas_int* /*lda*/,
      const number* /*vl*/,
      const number* /*vu*/,
      const types::blas_int* /*il*/,
      const types::blas_int* /*iu*/,
      const number* /*abstol*/,
      types::blas_int* /*m*/,
      number* /*w*/,
      number* /*z*/,
      const types::blas_int* /*ldz*/,
      types::blas_int* /*isuppz*/,
      number* /*work*/,
      types::blas_int* /*lwork*/,
      types::blas_int* /*iwork*/,
      types::blas_int* /*liwork*/,
      types::blas_int* /*info*/)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
syevr(const char*            jobz,
      const char*            range,
      const char*            uplo,
      const types::blas_int* n,
      double*                A,
      const types::blas_int* lda,
      const double*          vl,
      const double*          vu,
      const types::blas_int* il,
      const types::blas_int* iu,
      const double*          abstol,
      types::blas_int*       m,
      double*                w,
      double*                z,
      const types::blas_int* ldz,
      types::blas_int*       isuppz,
      double*                work,
      types::blas_int*       lwork,
      types::blas_int*       iwork,
      types::blas_int*       liwork,
      types::blas_int*       info)
{
  /*
   * Netlib and Atlas Lapack perform floating point tests (e.g. divide-by-zero) within the call to dsyevr
   * causing floating point exceptions to be thrown (at least in debug mode). Therefore, we wrap the calls
   * to dsyevr into the following code to suppress the exception.
   */
#  ifdef DEAL_II_HAVE_FP_EXCEPTIONS
  fenv_t fp_exceptions;
  feholdexcept(&fp_exceptions);
#  endif

  dsyevr_(jobz,
          range,
          uplo,
          n,
          A,
          lda,
          vl,
          vu,
          il,
          iu,
          abstol,
          m,
          w,
          z,
          ldz,
          isuppz,
          work,
          lwork,
          iwork,
          liwork,
          info);

#  ifdef DEAL_II_HAVE_FP_EXCEPTIONS
  fesetenv(&fp_exceptions);
#  endif
}

inline void
syevr(const char*            jobz,
      const char*            range,
      const char*            uplo,
      const types::blas_int* n,
      float*                 A,
      const types::blas_int* lda,
      const float*           vl,
      const float*           vu,
      const types::blas_int* il,
      const types::blas_int* iu,
      const float*           abstol,
      types::blas_int*       m,
      float*                 w,
      float*                 z,
      const types::blas_int* ldz,
      types::blas_int*       isuppz,
      float*                 work,
      types::blas_int*       lwork,
      types::blas_int*       iwork,
      types::blas_int*       liwork,
      types::blas_int*       info)
{
  /*
   * Netlib and Atlas Lapack perform floating point tests (e.g. divide-by-zero) within the call to ssyevr
   * causing floating point exceptions to be thrown (at least in debug mode). Therefore, we wrap the calls
   * to ssyevr into the following code to suppress the exception.
   */
#  ifdef DEAL_II_HAVE_FP_EXCEPTIONS
  fenv_t fp_exceptions;
  feholdexcept(&fp_exceptions);
#  endif

  ssyevr_(jobz,
          range,
          uplo,
          n,
          A,
          lda,
          vl,
          vu,
          il,
          iu,
          abstol,
          m,
          w,
          z,
          ldz,
          isuppz,
          work,
          lwork,
          iwork,
          liwork,
          info);

#  ifdef DEAL_II_HAVE_FP_EXCEPTIONS
  fesetenv(&fp_exceptions);
#  endif
}
#endif

/// Template wrapper for LAPACK functions dsygv and ssygv
template <typename number1,
          typename number2,
          typename number3,
          typename number4>
inline void
sygv(const types::blas_int*,
     const char*,
     const char*,
     const types::blas_int*,
     number1*,
     const types::blas_int*,
     number2*,
     const types::blas_int*,
     number3*,
     number4*,
     const types::blas_int*,
     types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
sygv(const types::blas_int* itype,
     const char*            jobz,
     const char*            uplo,
     const types::blas_int* n,
     double*                A,
     const types::blas_int* lda,
     double*                B,
     const types::blas_int* ldb,
     double*                w,
     double*                work,
     const types::blas_int* lwork,
     types::blas_int*       info)
{
  dsygv_(itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, info);
}
#else
inline void
sygv(const types::blas_int*,
     const char*,
     const char*,
     const types::blas_int*,
     double*,
     const types::blas_int*,
     double*,
     const types::blas_int*,
     double*,
     double*,
     const types::blas_int*,
     types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dsygv"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
sygv(const types::blas_int* itype,
     const char*            jobz,
     const char*            uplo,
     const types::blas_int* n,
     float*                 A,
     const types::blas_int* lda,
     float*                 B,
     const types::blas_int* ldb,
     float*                 w,
     float*                 work,
     const types::blas_int* lwork,
     types::blas_int*       info)
{
  ssygv_(itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, info);
}
#else
inline void
sygv(const types::blas_int*,
     const char*,
     const char*,
     const types::blas_int*,
     float*,
     const types::blas_int*,
     float*,
     const types::blas_int*,
     float*,
     float*,
     const types::blas_int*,
     types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("ssygv"));
}
#endif

/// Template wrapper for LAPACK functions dsygvx and ssygvx
template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5,
          typename number6,
          typename number7,
          typename number8>
inline void
sygvx(const types::blas_int*,
      const char*,
      const char*,
      const char*,
      const types::blas_int*,
      number1*,
      const types::blas_int*,
      number2*,
      const types::blas_int*,
      const number3*,
      const number4*,
      const types::blas_int*,
      const types::blas_int*,
      const number5*,
      types::blas_int*,
      number6*,
      number7*,
      const types::blas_int*,
      number8*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
sygvx(const types::blas_int* itype,
      const char*            jobz,
      const char*            range,
      const char*            uplo,
      const types::blas_int* n,
      double*                A,
      const types::blas_int* lda,
      double*                B,
      const types::blas_int* ldb,
      const double*          vl,
      const double*          vu,
      const types::blas_int* il,
      const types::blas_int* iu,
      const double*          abstol,
      types::blas_int*       m,
      double*                w,
      double*                z,
      const types::blas_int* ldz,
      double*                work,
      const types::blas_int* lwork,
      types::blas_int*       iwork,
      types::blas_int*       ifail,
      types::blas_int*       info)
{
  dsygvx_(itype,
          jobz,
          range,
          uplo,
          n,
          A,
          lda,
          B,
          ldb,
          vl,
          vu,
          il,
          iu,
          abstol,
          m,
          w,
          z,
          ldz,
          work,
          lwork,
          iwork,
          ifail,
          info);
}
#else
inline void
sygvx(const types::blas_int*,
      const char*,
      const char*,
      const char*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      const double*,
      const double*,
      const types::blas_int*,
      const types::blas_int*,
      const double*,
      types::blas_int*,
      double*,
      double*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dsygvx"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
sygvx(const types::blas_int* itype,
      const char*            jobz,
      const char*            range,
      const char*            uplo,
      const types::blas_int* n,
      float*                 A,
      const types::blas_int* lda,
      float*                 B,
      const types::blas_int* ldb,
      const float*           vl,
      const float*           vu,
      const types::blas_int* il,
      const types::blas_int* iu,
      const float*           abstol,
      types::blas_int*       m,
      float*                 w,
      float*                 z,
      const types::blas_int* ldz,
      float*                 work,
      const types::blas_int* lwork,
      types::blas_int*       iwork,
      types::blas_int*       ifail,
      types::blas_int*       info)
{
  ssygvx_(itype,
          jobz,
          range,
          uplo,
          n,
          A,
          lda,
          B,
          ldb,
          vl,
          vu,
          il,
          iu,
          abstol,
          m,
          w,
          z,
          ldz,
          work,
          lwork,
          iwork,
          ifail,
          info);
}
#else
inline void
sygvx(const types::blas_int*,
      const char*,
      const char*,
      const char*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      const float*,
      const float*,
      const types::blas_int*,
      const types::blas_int*,
      const float*,
      types::blas_int*,
      float*,
      float*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("ssygvx"));
}
#endif

/// Template wrapper for LAPACK functions dgesdd and sgesdd
template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5>
inline void
gesdd(const char*,
      const types::blas_int*,
      const types::blas_int*,
      number1*,
      const types::blas_int*,
      number2*,
      number3*,
      const types::blas_int*,
      number4*,
      const types::blas_int*,
      number5*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
gesdd(const char*            jobz,
      const types::blas_int* m,
      const types::blas_int* n,
      double*                A,
      const types::blas_int* lda,
      double*                s,
      double*                u,
      const types::blas_int* ldu,
      double*                vt,
      const types::blas_int* ldvt,
      double*                work,
      const types::blas_int* lwork,
      types::blas_int*       iwork,
      types::blas_int*       info)
{
  dgesdd_(jobz, m, n, A, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
}
#else
inline void
gesdd(const char*,
      const types::blas_int*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      double*,
      double*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dgesdd"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
gesdd(const char*            jobz,
      const types::blas_int* m,
      const types::blas_int* n,
      float*                 A,
      const types::blas_int* lda,
      float*                 s,
      float*                 u,
      const types::blas_int* ldu,
      float*                 vt,
      const types::blas_int* ldvt,
      float*                 work,
      const types::blas_int* lwork,
      types::blas_int*       iwork,
      types::blas_int*       info)
{
  sgesdd_(jobz, m, n, A, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
}
#else
inline void
gesdd(const char*,
      const types::blas_int*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      float*,
      float*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("sgesdd"));
}
#endif

/// Template wrapper for LAPACK functions dgesvd and sgesvd
template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5>
inline void
gesvd(types::blas_int*,
      types::blas_int*,
      const types::blas_int*,
      const types::blas_int*,
      number1*,
      const types::blas_int*,
      number2*,
      number3*,
      const types::blas_int*,
      number4*,
      const types::blas_int*,
      number5*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
gesvd(types::blas_int*       jobu,
      types::blas_int*       jobvt,
      const types::blas_int* n,
      const types::blas_int* m,
      double*                A,
      const types::blas_int* lda,
      double*                s,
      double*                u,
      const types::blas_int* ldu,
      double*                vt,
      const types::blas_int* ldvt,
      double*                work,
      const types::blas_int* lwork,
      types::blas_int*       info)
{
  dgesvd_(jobu, jobvt, n, m, A, lda, s, u, ldu, vt, ldvt, work, lwork, info);
}
#else
inline void
gesvd(types::blas_int*,
      types::blas_int*,
      const types::blas_int*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      double*,
      double*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dgesvd"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
gesvd(types::blas_int*       jobu,
      types::blas_int*       jobvt,
      const types::blas_int* n,
      const types::blas_int* m,
      float*                 A,
      const types::blas_int* lda,
      float*                 s,
      float*                 u,
      const types::blas_int* ldu,
      float*                 vt,
      const types::blas_int* ldvt,
      float*                 work,
      const types::blas_int* lwork,
      types::blas_int*       info)
{
  sgesvd_(jobu, jobvt, n, m, A, lda, s, u, ldu, vt, ldvt, work, lwork, info);
}
#else
inline void
gesvd(types::blas_int*,
      types::blas_int*,
      const types::blas_int*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      float*,
      float*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("sgesvd"));
}
#endif

/// Template wrapper for LAPACK functions dgelsd and sgelsd
template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5>
inline void
gelsd(const types::blas_int*,
      const types::blas_int*,
      const types::blas_int*,
      const number1*,
      const types::blas_int*,
      number2*,
      const types::blas_int*,
      number3*,
      const number4*,
      types::blas_int*,
      number5*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
gelsd(const types::blas_int* m,
      const types::blas_int* n,
      const types::blas_int* nrhs,
      const double*          A,
      const types::blas_int* lda,
      double*                B,
      const types::blas_int* ldb,
      double*                s,
      const double*          rcond,
      types::blas_int*       rank,
      double*                work,
      const types::blas_int* lwork,
      types::blas_int*       iwork,
      types::blas_int*       info)
{
  dgelsd_(m, n, nrhs, A, lda, B, ldb, s, rcond, rank, work, lwork, iwork, info);
}
#else
inline void
gelsd(const types::blas_int*,
      const types::blas_int*,
      const types::blas_int*,
      const double*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      double*,
      const double*,
      types::blas_int*,
      double*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dgelsd"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
gelsd(const types::blas_int* m,
      const types::blas_int* n,
      const types::blas_int* nrhs,
      const float*           A,
      const types::blas_int* lda,
      float*                 B,
      const types::blas_int* ldb,
      float*                 s,
      const float*           rcond,
      types::blas_int*       rank,
      float*                 work,
      const types::blas_int* lwork,
      types::blas_int*       iwork,
      types::blas_int*       info)
{
  sgelsd_(m, n, nrhs, A, lda, B, ldb, s, rcond, rank, work, lwork, iwork, info);
}
#else
inline void
gelsd(const types::blas_int*,
      const types::blas_int*,
      const types::blas_int*,
      const float*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      float*,
      const float*,
      types::blas_int*,
      float*,
      const types::blas_int*,
      types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("sgelsd"));
}
#endif

/// Template wrapper for LAPACK functions dstev and sstev
template <typename number1,
          typename number2,
          typename number3,
          typename number4>
inline void
stev(const char*,
     const types::blas_int*,
     number1*,
     number2*,
     number3*,
     const types::blas_int*,
     number4*,
     types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
stev(const char*            jobz,
     const types::blas_int* n,
     double*                d,
     double*                e,
     double*                z,
     const types::blas_int* ldz,
     double*                work,
     types::blas_int*       info)
{
  dstev_(jobz, n, d, e, z, ldz, work, info);
}
#else
inline void
stev(const char*,
     const types::blas_int*,
     double*,
     double*,
     double*,
     const types::blas_int*,
     double*,
     types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dstev"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
stev(const char*            jobz,
     const types::blas_int* n,
     float*                 d,
     float*                 e,
     float*                 z,
     const types::blas_int* ldz,
     float*                 work,
     types::blas_int*       info)
{
  sstev_(jobz, n, d, e, z, ldz, work, info);
}
#else
inline void
stev(const char*,
     const types::blas_int*,
     float*,
     float*,
     float*,
     const types::blas_int*,
     float*,
     types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("sstev"));
}
#endif

/// Template wrapper for LAPACK functions dlascl and slascl
template <typename number>
inline void
lascl(const char*,
      const types::blas_int*,
      const types::blas_int*,
      const number*,
      const number*,
      const types::blas_int*,
      const types::blas_int*,
      number*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
lascl(const char*            type,
      const types::blas_int* kl,
      const types::blas_int* ku,
      const double*          cfrom,
      const double*          cto,
      const types::blas_int* m,
      const types::blas_int* n,
      double*                A,
      const types::blas_int* lda,
      types::blas_int*       info)
{
  dlascl_(type, kl, ku, cfrom, cto, m, n, A, lda, info);
}
#else
inline void
lascl(const char*,
      const types::blas_int*,
      const types::blas_int*,
      const double*,
      const double*,
      const types::blas_int*,
      const types::blas_int*,
      double*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("dlascl"));
}
#endif

#ifdef DEAL_II_WITH_LAPACK
inline void
lascl(const char*            type,
      const types::blas_int* kl,
      const types::blas_int* ku,
      const float*           cfrom,
      const float*           cto,
      const types::blas_int* m,
      const types::blas_int* n,
      float*                 A,
      const types::blas_int* lda,
      types::blas_int*       info)
{
  slascl_(type, kl, ku, cfrom, cto, m, n, A, lda, info);
}
#else
inline void
lascl(const char*,
      const types::blas_int*,
      const types::blas_int*,
      const float*,
      const float*,
      const types::blas_int*,
      const types::blas_int*,
      float*,
      const types::blas_int*,
      types::blas_int*)
{
  Assert(false, LAPACKSupport::ExcMissing("slascl"));
}
#endif

/// Template wrapper for LAPACK functions dlamch and slamch
template <typename number>
inline void
lamch(const char* /*chmach*/, number& /*precision*/)
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_WITH_LAPACK
inline void
lamch(const char* chmach, double& precision)
{
  precision = dlamch_(chmach);
}

inline void
lamch(const char* chmach, float& precision)
{
  precision = slamch_(chmach);
}
#endif

DEAL_II_NAMESPACE_CLOSE

#endif
