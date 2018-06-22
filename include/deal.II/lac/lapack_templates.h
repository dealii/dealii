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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
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
  void
  saxpy_(const dealii::types::blas_int *n,
         const float *                  sa,
         const float *                  sx,
         const dealii::types::blas_int *incx,
         float *                        sy,
         const dealii::types::blas_int *incy);

  void
  daxpy_(const dealii::types::blas_int *n,
         const double *                 da,
         const double *                 dx,
         const dealii::types::blas_int *incx,
         double *                       dy,
         const dealii::types::blas_int *incy);

  void
  caxpy_(const dealii::types::blas_int *n,
         const std::complex<float> *    ca,
         const std::complex<float> *    cx,
         const dealii::types::blas_int *incx,
         std::complex<float> *          cy,
         const dealii::types::blas_int *incy);

  void
  zaxpy_(const dealii::types::blas_int *n,
         const std::complex<double> *   za,
         const std::complex<double> *   zx,
         const dealii::types::blas_int *incx,
         std::complex<double> *         zy,
         const dealii::types::blas_int *incy);

  void
  sgeev_(const char *                   jobvl,
         const char *                   jobvr,
         const dealii::types::blas_int *n,
         float *                        a,
         const dealii::types::blas_int *lda,
         float *                        wr,
         float *                        wi,
         float *                        vl,
         const dealii::types::blas_int *ldvl,
         float *                        vr,
         const dealii::types::blas_int *ldvr,
         float *                        work,
         const dealii::types::blas_int *lwork,
         dealii::types::blas_int *      info);

  void
  dgeev_(const char *                   jobvl,
         const char *                   jobvr,
         const dealii::types::blas_int *n,
         double *                       a,
         const dealii::types::blas_int *lda,
         double *                       wr,
         double *                       wi,
         double *                       vl,
         const dealii::types::blas_int *ldvl,
         double *                       vr,
         const dealii::types::blas_int *ldvr,
         double *                       work,
         const dealii::types::blas_int *lwork,
         dealii::types::blas_int *      info);

  void
  cgeev_(const char *                   jobvl,
         const char *                   jobvr,
         const dealii::types::blas_int *n,
         std::complex<float> *          a,
         const dealii::types::blas_int *lda,
         std::complex<float> *          w,
         std::complex<float> *          vl,
         const dealii::types::blas_int *ldvl,
         std::complex<float> *          vr,
         const dealii::types::blas_int *ldvr,
         std::complex<float> *          work,
         const dealii::types::blas_int *lwork,
         float *                        rwork,
         dealii::types::blas_int *      info);

  void
  zgeev_(const char *                   jobvl,
         const char *                   jobvr,
         const dealii::types::blas_int *n,
         std::complex<double> *         a,
         const dealii::types::blas_int *lda,
         std::complex<double> *         w,
         std::complex<double> *         vl,
         const dealii::types::blas_int *ldvl,
         std::complex<double> *         vr,
         const dealii::types::blas_int *ldvr,
         std::complex<double> *         work,
         const dealii::types::blas_int *lwork,
         double *                       rwork,
         dealii::types::blas_int *      info);

  void
  sgeevx_(const char *                   balanc,
          const char *                   jobvl,
          const char *                   jobvr,
          const char *                   sense,
          const dealii::types::blas_int *n,
          float *                        a,
          const dealii::types::blas_int *lda,
          float *                        wr,
          float *                        wi,
          float *                        vl,
          const dealii::types::blas_int *ldvl,
          float *                        vr,
          const dealii::types::blas_int *ldvr,
          dealii::types::blas_int *      ilo,
          dealii::types::blas_int *      ihi,
          float *                        scale,
          float *                        abnrm,
          float *                        rconde,
          float *                        rcondv,
          float *                        work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      info);

  void
  dgeevx_(const char *                   balanc,
          const char *                   jobvl,
          const char *                   jobvr,
          const char *                   sense,
          const dealii::types::blas_int *n,
          double *                       a,
          const dealii::types::blas_int *lda,
          double *                       wr,
          double *                       wi,
          double *                       vl,
          const dealii::types::blas_int *ldvl,
          double *                       vr,
          const dealii::types::blas_int *ldvr,
          dealii::types::blas_int *      ilo,
          dealii::types::blas_int *      ihi,
          double *                       scale,
          double *                       abnrm,
          double *                       rconde,
          double *                       rcondv,
          double *                       work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      info);

  void
  cgeevx_(const char *                   balanc,
          const char *                   jobvl,
          const char *                   jobvr,
          const char *                   sense,
          const dealii::types::blas_int *n,
          std::complex<float> *          a,
          const dealii::types::blas_int *lda,
          std::complex<float> *          w,
          std::complex<float> *          vl,
          const dealii::types::blas_int *ldvl,
          std::complex<float> *          vr,
          const dealii::types::blas_int *ldvr,
          dealii::types::blas_int *      ilo,
          dealii::types::blas_int *      ihi,
          float *                        scale,
          float *                        abnrm,
          float *                        rconde,
          float *                        rcondv,
          std::complex<float> *          work,
          const dealii::types::blas_int *lwork,
          float *                        rwork,
          dealii::types::blas_int *      info);

  void
  zgeevx_(const char *                   balanc,
          const char *                   jobvl,
          const char *                   jobvr,
          const char *                   sense,
          const dealii::types::blas_int *n,
          std::complex<double> *         a,
          const dealii::types::blas_int *lda,
          std::complex<double> *         w,
          std::complex<double> *         vl,
          const dealii::types::blas_int *ldvl,
          std::complex<double> *         vr,
          const dealii::types::blas_int *ldvr,
          dealii::types::blas_int *      ilo,
          dealii::types::blas_int *      ihi,
          double *                       scale,
          double *                       abnrm,
          double *                       rconde,
          double *                       rcondv,
          std::complex<double> *         work,
          const dealii::types::blas_int *lwork,
          double *                       rwork,
          dealii::types::blas_int *      info);

  void
  sgelsd_(const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          float *                        a,
          const dealii::types::blas_int *lda,
          float *                        b,
          const dealii::types::blas_int *ldb,
          float *                        s,
          const float *                  rcond,
          dealii::types::blas_int *      rank,
          float *                        work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      info);

  void
  dgelsd_(const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          double *                       a,
          const dealii::types::blas_int *lda,
          double *                       b,
          const dealii::types::blas_int *ldb,
          double *                       s,
          const double *                 rcond,
          dealii::types::blas_int *      rank,
          double *                       work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      info);

  void
  cgelsd_(const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          std::complex<float> *          a,
          const dealii::types::blas_int *lda,
          std::complex<float> *          b,
          const dealii::types::blas_int *ldb,
          float *                        s,
          const float *                  rcond,
          dealii::types::blas_int *      rank,
          std::complex<float> *          work,
          const dealii::types::blas_int *lwork,
          float *                        rwork,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      info);

  void
  zgelsd_(const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          std::complex<double> *         a,
          const dealii::types::blas_int *lda,
          std::complex<double> *         b,
          const dealii::types::blas_int *ldb,
          double *                       s,
          const double *                 rcond,
          dealii::types::blas_int *      rank,
          std::complex<double> *         work,
          const dealii::types::blas_int *lwork,
          double *                       rwork,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      info);

  void
  sgemm_(const char *                   transa,
         const char *                   transb,
         const dealii::types::blas_int *m,
         const dealii::types::blas_int *n,
         const dealii::types::blas_int *k,
         const float *                  alpha,
         const float *                  a,
         const dealii::types::blas_int *lda,
         const float *                  b,
         const dealii::types::blas_int *ldb,
         const float *                  beta,
         float *                        c,
         const dealii::types::blas_int *ldc);

  void
  dgemm_(const char *                   transa,
         const char *                   transb,
         const dealii::types::blas_int *m,
         const dealii::types::blas_int *n,
         const dealii::types::blas_int *k,
         const double *                 alpha,
         const double *                 a,
         const dealii::types::blas_int *lda,
         const double *                 b,
         const dealii::types::blas_int *ldb,
         const double *                 beta,
         double *                       c,
         const dealii::types::blas_int *ldc);

  void
  cgemm_(const char *                   transa,
         const char *                   transb,
         const dealii::types::blas_int *m,
         const dealii::types::blas_int *n,
         const dealii::types::blas_int *k,
         const std::complex<float> *    alpha,
         const std::complex<float> *    a,
         const dealii::types::blas_int *lda,
         const std::complex<float> *    b,
         const dealii::types::blas_int *ldb,
         const std::complex<float> *    beta,
         std::complex<float> *          c,
         const dealii::types::blas_int *ldc);

  void
  zgemm_(const char *                   transa,
         const char *                   transb,
         const dealii::types::blas_int *m,
         const dealii::types::blas_int *n,
         const dealii::types::blas_int *k,
         const std::complex<double> *   alpha,
         const std::complex<double> *   a,
         const dealii::types::blas_int *lda,
         const std::complex<double> *   b,
         const dealii::types::blas_int *ldb,
         const std::complex<double> *   beta,
         std::complex<double> *         c,
         const dealii::types::blas_int *ldc);

  void
  sgemv_(const char *                   trans,
         const dealii::types::blas_int *m,
         const dealii::types::blas_int *n,
         const float *                  alpha,
         const float *                  a,
         const dealii::types::blas_int *lda,
         const float *                  x,
         const dealii::types::blas_int *incx,
         const float *                  beta,
         float *                        y,
         const dealii::types::blas_int *incy);

  void
  dgemv_(const char *                   trans,
         const dealii::types::blas_int *m,
         const dealii::types::blas_int *n,
         const double *                 alpha,
         const double *                 a,
         const dealii::types::blas_int *lda,
         const double *                 x,
         const dealii::types::blas_int *incx,
         const double *                 beta,
         double *                       y,
         const dealii::types::blas_int *incy);

  void
  cgemv_(const char *                   trans,
         const dealii::types::blas_int *m,
         const dealii::types::blas_int *n,
         const std::complex<float> *    alpha,
         const std::complex<float> *    a,
         const dealii::types::blas_int *lda,
         const std::complex<float> *    x,
         const dealii::types::blas_int *incx,
         const std::complex<float> *    beta,
         std::complex<float> *          y,
         const dealii::types::blas_int *incy);

  void
  zgemv_(const char *                   trans,
         const dealii::types::blas_int *m,
         const dealii::types::blas_int *n,
         const std::complex<double> *   alpha,
         const std::complex<double> *   a,
         const dealii::types::blas_int *lda,
         const std::complex<double> *   x,
         const dealii::types::blas_int *incx,
         const std::complex<double> *   beta,
         std::complex<double> *         y,
         const dealii::types::blas_int *incy);

  void
  sgeqrf_(const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          float *                        a,
          const dealii::types::blas_int *lda,
          float *                        tau,
          float *                        work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      info);

  void
  dgeqrf_(const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          double *                       a,
          const dealii::types::blas_int *lda,
          double *                       tau,
          double *                       work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      info);

  void
  cgeqrf_(const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          std::complex<float> *          a,
          const dealii::types::blas_int *lda,
          std::complex<float> *          tau,
          std::complex<float> *          work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      info);

  void
  zgeqrf_(const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          std::complex<double> *         a,
          const dealii::types::blas_int *lda,
          std::complex<double> *         tau,
          std::complex<double> *         work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      info);

  void
  sgesdd_(const char *                   jobz,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          float *                        a,
          const dealii::types::blas_int *lda,
          float *                        s,
          float *                        u,
          const dealii::types::blas_int *ldu,
          float *                        vt,
          const dealii::types::blas_int *ldvt,
          float *                        work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      info);

  void
  dgesdd_(const char *                   jobz,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          double *                       a,
          const dealii::types::blas_int *lda,
          double *                       s,
          double *                       u,
          const dealii::types::blas_int *ldu,
          double *                       vt,
          const dealii::types::blas_int *ldvt,
          double *                       work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      info);

  void
  cgesdd_(const char *                   jobz,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          std::complex<float> *          a,
          const dealii::types::blas_int *lda,
          float *                        s,
          std::complex<float> *          u,
          const dealii::types::blas_int *ldu,
          std::complex<float> *          vt,
          const dealii::types::blas_int *ldvt,
          std::complex<float> *          work,
          const dealii::types::blas_int *lwork,
          float *                        rwork,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      info);

  void
  zgesdd_(const char *                   jobz,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          std::complex<double> *         a,
          const dealii::types::blas_int *lda,
          double *                       s,
          std::complex<double> *         u,
          const dealii::types::blas_int *ldu,
          std::complex<double> *         vt,
          const dealii::types::blas_int *ldvt,
          std::complex<double> *         work,
          const dealii::types::blas_int *lwork,
          double *                       rwork,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      info);

  void
  sgesvd_(const char *                   jobu,
          const char *                   jobvt,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          float *                        a,
          const dealii::types::blas_int *lda,
          float *                        s,
          float *                        u,
          const dealii::types::blas_int *ldu,
          float *                        vt,
          const dealii::types::blas_int *ldvt,
          float *                        work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      info);

  void
  dgesvd_(const char *                   jobu,
          const char *                   jobvt,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          double *                       a,
          const dealii::types::blas_int *lda,
          double *                       s,
          double *                       u,
          const dealii::types::blas_int *ldu,
          double *                       vt,
          const dealii::types::blas_int *ldvt,
          double *                       work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      info);

  void
  cgesvd_(const char *                   jobu,
          const char *                   jobvt,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          std::complex<float> *          a,
          const dealii::types::blas_int *lda,
          float *                        s,
          std::complex<float> *          u,
          const dealii::types::blas_int *ldu,
          std::complex<float> *          vt,
          const dealii::types::blas_int *ldvt,
          std::complex<float> *          work,
          const dealii::types::blas_int *lwork,
          float *                        rwork,
          dealii::types::blas_int *      info);

  void
  zgesvd_(const char *                   jobu,
          const char *                   jobvt,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          std::complex<double> *         a,
          const dealii::types::blas_int *lda,
          double *                       s,
          std::complex<double> *         u,
          const dealii::types::blas_int *ldu,
          std::complex<double> *         vt,
          const dealii::types::blas_int *ldvt,
          std::complex<double> *         work,
          const dealii::types::blas_int *lwork,
          double *                       rwork,
          dealii::types::blas_int *      info);

  void
  sgetrf_(const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          float *                        a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      ipiv,
          dealii::types::blas_int *      info);

  void
  dgetrf_(const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          double *                       a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      ipiv,
          dealii::types::blas_int *      info);

  void
  cgetrf_(const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          std::complex<float> *          a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      ipiv,
          dealii::types::blas_int *      info);

  void
  zgetrf_(const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          std::complex<double> *         a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      ipiv,
          dealii::types::blas_int *      info);

  void
  sgetri_(const dealii::types::blas_int *n,
          float *                        a,
          const dealii::types::blas_int *lda,
          const dealii::types::blas_int *ipiv,
          float *                        work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      info);

  void
  dgetri_(const dealii::types::blas_int *n,
          double *                       a,
          const dealii::types::blas_int *lda,
          const dealii::types::blas_int *ipiv,
          double *                       work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      info);

  void
  cgetri_(const dealii::types::blas_int *n,
          std::complex<float> *          a,
          const dealii::types::blas_int *lda,
          const dealii::types::blas_int *ipiv,
          std::complex<float> *          work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      info);

  void
  zgetri_(const dealii::types::blas_int *n,
          std::complex<double> *         a,
          const dealii::types::blas_int *lda,
          const dealii::types::blas_int *ipiv,
          std::complex<double> *         work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      info);

  void
  sgetrs_(const char *                   trans,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          const float *                  a,
          const dealii::types::blas_int *lda,
          const dealii::types::blas_int *ipiv,
          float *                        b,
          const dealii::types::blas_int *ldb,
          dealii::types::blas_int *      info);

  void
  dgetrs_(const char *                   trans,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          const double *                 a,
          const dealii::types::blas_int *lda,
          const dealii::types::blas_int *ipiv,
          double *                       b,
          const dealii::types::blas_int *ldb,
          dealii::types::blas_int *      info);

  void
  cgetrs_(const char *                   trans,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          const std::complex<float> *    a,
          const dealii::types::blas_int *lda,
          const dealii::types::blas_int *ipiv,
          std::complex<float> *          b,
          const dealii::types::blas_int *ldb,
          dealii::types::blas_int *      info);

  void
  zgetrs_(const char *                   trans,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          const std::complex<double> *   a,
          const dealii::types::blas_int *lda,
          const dealii::types::blas_int *ipiv,
          std::complex<double> *         b,
          const dealii::types::blas_int *ldb,
          dealii::types::blas_int *      info);

  float
  slamch_(const char *cmach);

  double
  dlamch_(const char *cmach);

  float
  slange_(const char *                   norm,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          const float *                  a,
          const dealii::types::blas_int *lda,
          float *                        work);

  double
  dlange_(const char *                   norm,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          const double *                 a,
          const dealii::types::blas_int *lda,
          double *                       work);

  float
  clange_(const char *                   norm,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          const std::complex<float> *    a,
          const dealii::types::blas_int *lda,
          float *                        work);

  double
  zlange_(const char *                   norm,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          const std::complex<double> *   a,
          const dealii::types::blas_int *lda,
          double *                       work);

  float
  slansy_(const char *                   norm,
          const char *                   uplo,
          const dealii::types::blas_int *n,
          const float *                  a,
          const dealii::types::blas_int *lda,
          float *                        work);

  double
  dlansy_(const char *                   norm,
          const char *                   uplo,
          const dealii::types::blas_int *n,
          const double *                 a,
          const dealii::types::blas_int *lda,
          double *                       work);

  float
  clansy_(const char *                   norm,
          const char *                   uplo,
          const dealii::types::blas_int *n,
          const std::complex<float> *    a,
          const dealii::types::blas_int *lda,
          float *                        work);

  double
  zlansy_(const char *                   norm,
          const char *                   uplo,
          const dealii::types::blas_int *n,
          const std::complex<double> *   a,
          const dealii::types::blas_int *lda,
          double *                       work);

  void
  slascl_(const char *                   type,
          const dealii::types::blas_int *kl,
          const dealii::types::blas_int *ku,
          const float *                  cfrom,
          const float *                  cto,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          float *                        a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      info);

  void
  dlascl_(const char *                   type,
          const dealii::types::blas_int *kl,
          const dealii::types::blas_int *ku,
          const double *                 cfrom,
          const double *                 cto,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          double *                       a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      info);

  void
  clascl_(const char *                   type,
          const dealii::types::blas_int *kl,
          const dealii::types::blas_int *ku,
          const float *                  cfrom,
          const float *                  cto,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          std::complex<float> *          a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      info);

  void
  zlascl_(const char *                   type,
          const dealii::types::blas_int *kl,
          const dealii::types::blas_int *ku,
          const double *                 cfrom,
          const double *                 cto,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          std::complex<double> *         a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      info);

  void
  sorgqr_(const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *k,
          float *                        a,
          const dealii::types::blas_int *lda,
          const float *                  tau,
          float *                        work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      info);

  void
  dorgqr_(const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *k,
          double *                       a,
          const dealii::types::blas_int *lda,
          const double *                 tau,
          double *                       work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      info);

  void
  sormqr_(const char *                   side,
          const char *                   trans,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *k,
          const float *                  a,
          const dealii::types::blas_int *lda,
          const float *                  tau,
          float *                        c,
          const dealii::types::blas_int *ldc,
          float *                        work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      info);

  void
  dormqr_(const char *                   side,
          const char *                   trans,
          const dealii::types::blas_int *m,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *k,
          const double *                 a,
          const dealii::types::blas_int *lda,
          const double *                 tau,
          double *                       c,
          const dealii::types::blas_int *ldc,
          double *                       work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      info);

  void
  spocon_(const char *                   uplo,
          const dealii::types::blas_int *n,
          const float *                  a,
          const dealii::types::blas_int *lda,
          const float *                  anorm,
          float *                        rcond,
          float *                        work,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      info);

  void
  dpocon_(const char *                   uplo,
          const dealii::types::blas_int *n,
          const double *                 a,
          const dealii::types::blas_int *lda,
          const double *                 anorm,
          double *                       rcond,
          double *                       work,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      info);

  void
  cpocon_(const char *                   uplo,
          const dealii::types::blas_int *n,
          const std::complex<float> *    a,
          const dealii::types::blas_int *lda,
          const float *                  anorm,
          float *                        rcond,
          std::complex<float> *          work,
          float *                        rwork,
          dealii::types::blas_int *      info);

  void
  zpocon_(const char *                   uplo,
          const dealii::types::blas_int *n,
          const std::complex<double> *   a,
          const dealii::types::blas_int *lda,
          const double *                 anorm,
          double *                       rcond,
          std::complex<double> *         work,
          double *                       rwork,
          dealii::types::blas_int *      info);

  void
  spotrf_(const char *                   uplo,
          const dealii::types::blas_int *n,
          float *                        a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      info);

  void
  dpotrf_(const char *                   uplo,
          const dealii::types::blas_int *n,
          double *                       a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      info);

  void
  cpotrf_(const char *                   uplo,
          const dealii::types::blas_int *n,
          std::complex<float> *          a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      info);

  void
  zpotrf_(const char *                   uplo,
          const dealii::types::blas_int *n,
          std::complex<double> *         a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      info);

  void
  spotri_(const char *                   uplo,
          const dealii::types::blas_int *n,
          float *                        a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      info);

  void
  dpotri_(const char *                   uplo,
          const dealii::types::blas_int *n,
          double *                       a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      info);

  void
  cpotri_(const char *                   uplo,
          const dealii::types::blas_int *n,
          std::complex<float> *          a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      info);

  void
  zpotri_(const char *                   uplo,
          const dealii::types::blas_int *n,
          std::complex<double> *         a,
          const dealii::types::blas_int *lda,
          dealii::types::blas_int *      info);

  void
  spotrs_(const char *                   uplo,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          const float *                  a,
          const dealii::types::blas_int *lda,
          float *                        b,
          const dealii::types::blas_int *ldb,
          dealii::types::blas_int *      info);

  void
  dpotrs_(const char *                   uplo,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          const double *                 a,
          const dealii::types::blas_int *lda,
          double *                       b,
          const dealii::types::blas_int *ldb,
          dealii::types::blas_int *      info);

  void
  cpotrs_(const char *                   uplo,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          const std::complex<float> *    a,
          const dealii::types::blas_int *lda,
          std::complex<float> *          b,
          const dealii::types::blas_int *ldb,
          dealii::types::blas_int *      info);

  void
  zpotrs_(const char *                   uplo,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          const std::complex<double> *   a,
          const dealii::types::blas_int *lda,
          std::complex<double> *         b,
          const dealii::types::blas_int *ldb,
          dealii::types::blas_int *      info);

  void
  sstev_(const char *                   jobz,
         const dealii::types::blas_int *n,
         float *                        d,
         float *                        e,
         float *                        z,
         const dealii::types::blas_int *ldz,
         float *                        work,
         dealii::types::blas_int *      info);

  void
  dstev_(const char *                   jobz,
         const dealii::types::blas_int *n,
         double *                       d,
         double *                       e,
         double *                       z,
         const dealii::types::blas_int *ldz,
         double *                       work,
         dealii::types::blas_int *      info);

  void
  ssyev_(const char *                   jobz,
         const char *                   uplo,
         const dealii::types::blas_int *n,
         float *                        a,
         const dealii::types::blas_int *lda,
         float *                        w,
         float *                        work,
         const dealii::types::blas_int *lwork,
         dealii::types::blas_int *      info);

  void
  dsyev_(const char *                   jobz,
         const char *                   uplo,
         const dealii::types::blas_int *n,
         double *                       a,
         const dealii::types::blas_int *lda,
         double *                       w,
         double *                       work,
         const dealii::types::blas_int *lwork,
         dealii::types::blas_int *      info);

  void
  ssyevr_(const char *                   jobz,
          const char *                   range,
          const char *                   uplo,
          const dealii::types::blas_int *n,
          float *                        a,
          const dealii::types::blas_int *lda,
          const float *                  vl,
          const float *                  vu,
          const dealii::types::blas_int *il,
          const dealii::types::blas_int *iu,
          const float *                  abstol,
          dealii::types::blas_int *      m,
          float *                        w,
          float *                        z,
          const dealii::types::blas_int *ldz,
          dealii::types::blas_int *      isuppz,
          float *                        work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      iwork,
          const dealii::types::blas_int *liwork,
          dealii::types::blas_int *      info);

  void
  dsyevr_(const char *                   jobz,
          const char *                   range,
          const char *                   uplo,
          const dealii::types::blas_int *n,
          double *                       a,
          const dealii::types::blas_int *lda,
          const double *                 vl,
          const double *                 vu,
          const dealii::types::blas_int *il,
          const dealii::types::blas_int *iu,
          const double *                 abstol,
          dealii::types::blas_int *      m,
          double *                       w,
          double *                       z,
          const dealii::types::blas_int *ldz,
          dealii::types::blas_int *      isuppz,
          double *                       work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      iwork,
          const dealii::types::blas_int *liwork,
          dealii::types::blas_int *      info);

  void
  ssyevx_(const char *                   jobz,
          const char *                   range,
          const char *                   uplo,
          const dealii::types::blas_int *n,
          float *                        a,
          const dealii::types::blas_int *lda,
          const float *                  vl,
          const float *                  vu,
          const dealii::types::blas_int *il,
          const dealii::types::blas_int *iu,
          const float *                  abstol,
          dealii::types::blas_int *      m,
          float *                        w,
          float *                        z,
          const dealii::types::blas_int *ldz,
          float *                        work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      ifail,
          dealii::types::blas_int *      info);

  void
  dsyevx_(const char *                   jobz,
          const char *                   range,
          const char *                   uplo,
          const dealii::types::blas_int *n,
          double *                       a,
          const dealii::types::blas_int *lda,
          const double *                 vl,
          const double *                 vu,
          const dealii::types::blas_int *il,
          const dealii::types::blas_int *iu,
          const double *                 abstol,
          dealii::types::blas_int *      m,
          double *                       w,
          double *                       z,
          const dealii::types::blas_int *ldz,
          double *                       work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      ifail,
          dealii::types::blas_int *      info);

  void
  ssygv_(const dealii::types::blas_int *itype,
         const char *                   jobz,
         const char *                   uplo,
         const dealii::types::blas_int *n,
         float *                        a,
         const dealii::types::blas_int *lda,
         float *                        b,
         const dealii::types::blas_int *ldb,
         float *                        w,
         float *                        work,
         const dealii::types::blas_int *lwork,
         dealii::types::blas_int *      info);

  void
  dsygv_(const dealii::types::blas_int *itype,
         const char *                   jobz,
         const char *                   uplo,
         const dealii::types::blas_int *n,
         double *                       a,
         const dealii::types::blas_int *lda,
         double *                       b,
         const dealii::types::blas_int *ldb,
         double *                       w,
         double *                       work,
         const dealii::types::blas_int *lwork,
         dealii::types::blas_int *      info);

  void
  ssygvx_(const dealii::types::blas_int *itype,
          const char *                   jobz,
          const char *                   range,
          const char *                   uplo,
          const dealii::types::blas_int *n,
          float *                        a,
          const dealii::types::blas_int *lda,
          float *                        b,
          const dealii::types::blas_int *ldb,
          const float *                  vl,
          const float *                  vu,
          const dealii::types::blas_int *il,
          const dealii::types::blas_int *iu,
          const float *                  abstol,
          dealii::types::blas_int *      m,
          float *                        w,
          float *                        z,
          const dealii::types::blas_int *ldz,
          float *                        work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      ifail,
          dealii::types::blas_int *      info);

  void
  dsygvx_(const dealii::types::blas_int *itype,
          const char *                   jobz,
          const char *                   range,
          const char *                   uplo,
          const dealii::types::blas_int *n,
          double *                       a,
          const dealii::types::blas_int *lda,
          double *                       b,
          const dealii::types::blas_int *ldb,
          const double *                 vl,
          const double *                 vu,
          const dealii::types::blas_int *il,
          const dealii::types::blas_int *iu,
          const double *                 abstol,
          dealii::types::blas_int *      m,
          double *                       w,
          double *                       z,
          const dealii::types::blas_int *ldz,
          double *                       work,
          const dealii::types::blas_int *lwork,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      ifail,
          dealii::types::blas_int *      info);

  void
  ssyr_(const char *                   uplo,
        const dealii::types::blas_int *n,
        const float *                  alpha,
        const float *                  x,
        const dealii::types::blas_int *incx,
        float *                        a,
        const dealii::types::blas_int *lda);

  void
  dsyr_(const char *                   uplo,
        const dealii::types::blas_int *n,
        const double *                 alpha,
        const double *                 x,
        const dealii::types::blas_int *incx,
        double *                       a,
        const dealii::types::blas_int *lda);

  void
  ssyrk_(const char *                   uplo,
         const char *                   trans,
         const dealii::types::blas_int *n,
         const dealii::types::blas_int *k,
         const float *                  alpha,
         const float *                  a,
         const dealii::types::blas_int *lda,
         const float *                  beta,
         float *                        c,
         const dealii::types::blas_int *ldc);

  void
  dsyrk_(const char *                   uplo,
         const char *                   trans,
         const dealii::types::blas_int *n,
         const dealii::types::blas_int *k,
         const double *                 alpha,
         const double *                 a,
         const dealii::types::blas_int *lda,
         const double *                 beta,
         double *                       c,
         const dealii::types::blas_int *ldc);

  void
  csyrk_(const char *                   uplo,
         const char *                   trans,
         const dealii::types::blas_int *n,
         const dealii::types::blas_int *k,
         const std::complex<float> *    alpha,
         const std::complex<float> *    a,
         const dealii::types::blas_int *lda,
         const std::complex<float> *    beta,
         std::complex<float> *          c,
         const dealii::types::blas_int *ldc);

  void
  zsyrk_(const char *                   uplo,
         const char *                   trans,
         const dealii::types::blas_int *n,
         const dealii::types::blas_int *k,
         const std::complex<double> *   alpha,
         const std::complex<double> *   a,
         const dealii::types::blas_int *lda,
         const std::complex<double> *   beta,
         std::complex<double> *         c,
         const dealii::types::blas_int *ldc);

  void
  strcon_(const char *                   norm,
          const char *                   uplo,
          const char *                   diag,
          const dealii::types::blas_int *n,
          const float *                  a,
          const dealii::types::blas_int *lda,
          float *                        rcond,
          float *                        work,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      info);

  void
  dtrcon_(const char *                   norm,
          const char *                   uplo,
          const char *                   diag,
          const dealii::types::blas_int *n,
          const double *                 a,
          const dealii::types::blas_int *lda,
          double *                       rcond,
          double *                       work,
          dealii::types::blas_int *      iwork,
          dealii::types::blas_int *      info);

  void
  ctrcon_(const char *                   norm,
          const char *                   uplo,
          const char *                   diag,
          const dealii::types::blas_int *n,
          const std::complex<float> *    a,
          const dealii::types::blas_int *lda,
          float *                        rcond,
          std::complex<float> *          work,
          float *                        rwork,
          dealii::types::blas_int *      info);

  void
  ztrcon_(const char *                   norm,
          const char *                   uplo,
          const char *                   diag,
          const dealii::types::blas_int *n,
          const std::complex<double> *   a,
          const dealii::types::blas_int *lda,
          double *                       rcond,
          std::complex<double> *         work,
          double *                       rwork,
          dealii::types::blas_int *      info);

  void
  strmv_(const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const dealii::types::blas_int *n,
         const float *                  a,
         const dealii::types::blas_int *lda,
         float *                        x,
         const dealii::types::blas_int *incx);

  void
  dtrmv_(const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const dealii::types::blas_int *n,
         const double *                 a,
         const dealii::types::blas_int *lda,
         double *                       x,
         const dealii::types::blas_int *incx);

  void
  ctrmv_(const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const dealii::types::blas_int *n,
         const std::complex<float> *    a,
         const dealii::types::blas_int *lda,
         std::complex<float> *          x,
         const dealii::types::blas_int *incx);

  void
  ztrmv_(const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const dealii::types::blas_int *n,
         const std::complex<double> *   a,
         const dealii::types::blas_int *lda,
         std::complex<double> *         x,
         const dealii::types::blas_int *incx);

  void
  strtrs_(const char *                   uplo,
          const char *                   trans,
          const char *                   diag,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          const float *                  a,
          const dealii::types::blas_int *lda,
          float *                        b,
          const dealii::types::blas_int *ldb,
          dealii::types::blas_int *      info);

  void
  dtrtrs_(const char *                   uplo,
          const char *                   trans,
          const char *                   diag,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          const double *                 a,
          const dealii::types::blas_int *lda,
          double *                       b,
          const dealii::types::blas_int *ldb,
          dealii::types::blas_int *      info);

  void
  ctrtrs_(const char *                   uplo,
          const char *                   trans,
          const char *                   diag,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          const std::complex<float> *    a,
          const dealii::types::blas_int *lda,
          std::complex<float> *          b,
          const dealii::types::blas_int *ldb,
          dealii::types::blas_int *      info);

  void
  ztrtrs_(const char *                   uplo,
          const char *                   trans,
          const char *                   diag,
          const dealii::types::blas_int *n,
          const dealii::types::blas_int *nrhs,
          const std::complex<double> *   a,
          const dealii::types::blas_int *lda,
          std::complex<double> *         b,
          const dealii::types::blas_int *ldb,
          dealii::types::blas_int *      info);
}

DEAL_II_NAMESPACE_OPEN



template <typename number1, typename number2, typename number3>
inline void
axpy(const dealii::types::blas_int *,
     const number1 *,
     const number2 *,
     const dealii::types::blas_int *,
     number3 *,
     const dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
axpy(const dealii::types::blas_int *n,
     const float *                  sa,
     const float *                  sx,
     const dealii::types::blas_int *incx,
     float *                        sy,
     const dealii::types::blas_int *incy)
{
#ifdef DEAL_II_WITH_LAPACK
  saxpy_(n, sa, sx, incx, sy, incy);
#else
  (void)n;
  (void)sa;
  (void)sx;
  (void)incx;
  (void)sy;
  (void)incy;
  Assert(false, LAPACKSupport::ExcMissing("saxpy"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
axpy(const dealii::types::blas_int *n,
     const double *                 da,
     const double *                 dx,
     const dealii::types::blas_int *incx,
     double *                       dy,
     const dealii::types::blas_int *incy)
{
#ifdef DEAL_II_WITH_LAPACK
  daxpy_(n, da, dx, incx, dy, incy);
#else
  (void)n;
  (void)da;
  (void)dx;
  (void)incx;
  (void)dy;
  (void)incy;
  Assert(false, LAPACKSupport::ExcMissing("daxpy"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
axpy(const dealii::types::blas_int *n,
     const std::complex<float> *    ca,
     const std::complex<float> *    cx,
     const dealii::types::blas_int *incx,
     std::complex<float> *          cy,
     const dealii::types::blas_int *incy)
{
#ifdef DEAL_II_WITH_LAPACK
  caxpy_(n, ca, cx, incx, cy, incy);
#else
  (void)n;
  (void)ca;
  (void)cx;
  (void)incx;
  (void)cy;
  (void)incy;
  Assert(false, LAPACKSupport::ExcMissing("caxpy"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
axpy(const dealii::types::blas_int *n,
     const std::complex<double> *   za,
     const std::complex<double> *   zx,
     const dealii::types::blas_int *incx,
     std::complex<double> *         zy,
     const dealii::types::blas_int *incy)
{
#ifdef DEAL_II_WITH_LAPACK
  zaxpy_(n, za, zx, incx, zy, incy);
#else
  (void)n;
  (void)za;
  (void)zx;
  (void)incx;
  (void)zy;
  (void)incy;
  Assert(false, LAPACKSupport::ExcMissing("zaxpy"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5,
          typename number6>
inline void
geev(const char *,
     const char *,
     const dealii::types::blas_int *,
     number1 *,
     const dealii::types::blas_int *,
     number2 *,
     number3 *,
     number4 *,
     const dealii::types::blas_int *,
     number5 *,
     const dealii::types::blas_int *,
     number6 *,
     const dealii::types::blas_int *,
     dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
geev(const char *                   jobvl,
     const char *                   jobvr,
     const dealii::types::blas_int *n,
     float *                        a,
     const dealii::types::blas_int *lda,
     float *                        wr,
     float *                        wi,
     float *                        vl,
     const dealii::types::blas_int *ldvl,
     float *                        vr,
     const dealii::types::blas_int *ldvr,
     float *                        work,
     const dealii::types::blas_int *lwork,
     dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  sgeev_(
    jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
#else
  (void)jobvl;
  (void)jobvr;
  (void)n;
  (void)a;
  (void)lda;
  (void)wr;
  (void)wi;
  (void)vl;
  (void)ldvl;
  (void)vr;
  (void)ldvr;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("sgeev"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
geev(const char *                   jobvl,
     const char *                   jobvr,
     const dealii::types::blas_int *n,
     double *                       a,
     const dealii::types::blas_int *lda,
     double *                       wr,
     double *                       wi,
     double *                       vl,
     const dealii::types::blas_int *ldvl,
     double *                       vr,
     const dealii::types::blas_int *ldvr,
     double *                       work,
     const dealii::types::blas_int *lwork,
     dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dgeev_(
    jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
#else
  (void)jobvl;
  (void)jobvr;
  (void)n;
  (void)a;
  (void)lda;
  (void)wr;
  (void)wi;
  (void)vl;
  (void)ldvl;
  (void)vr;
  (void)ldvr;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dgeev"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
geev(const char *                   jobvl,
     const char *                   jobvr,
     const dealii::types::blas_int *n,
     std::complex<float> *          a,
     const dealii::types::blas_int *lda,
     std::complex<float> *          w,
     std::complex<float> *          vl,
     const dealii::types::blas_int *ldvl,
     std::complex<float> *          vr,
     const dealii::types::blas_int *ldvr,
     std::complex<float> *          work,
     const dealii::types::blas_int *lwork,
     float *                        rwork,
     dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  cgeev_(
    jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
#else
  (void)jobvl;
  (void)jobvr;
  (void)n;
  (void)a;
  (void)lda;
  (void)w;
  (void)vl;
  (void)ldvl;
  (void)vr;
  (void)ldvr;
  (void)work;
  (void)lwork;
  (void)rwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("cgeev"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
geev(const char *                   jobvl,
     const char *                   jobvr,
     const dealii::types::blas_int *n,
     std::complex<double> *         a,
     const dealii::types::blas_int *lda,
     std::complex<double> *         w,
     std::complex<double> *         vl,
     const dealii::types::blas_int *ldvl,
     std::complex<double> *         vr,
     const dealii::types::blas_int *ldvr,
     std::complex<double> *         work,
     const dealii::types::blas_int *lwork,
     double *                       rwork,
     dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  zgeev_(
    jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
#else
  (void)jobvl;
  (void)jobvr;
  (void)n;
  (void)a;
  (void)lda;
  (void)w;
  (void)vl;
  (void)ldvl;
  (void)vr;
  (void)ldvr;
  (void)work;
  (void)lwork;
  (void)rwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("zgeev"));
#endif // DEAL_II_WITH_LAPACK
}



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
geevx(const char *,
      const char *,
      const char *,
      const char *,
      const dealii::types::blas_int *,
      number1 *,
      const dealii::types::blas_int *,
      number2 *,
      number3 *,
      number4 *,
      const dealii::types::blas_int *,
      number5 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *,
      dealii::types::blas_int *,
      number6 *,
      number7 *,
      number8 *,
      number9 *,
      number10 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
geevx(const char *                   balanc,
      const char *                   jobvl,
      const char *                   jobvr,
      const char *                   sense,
      const dealii::types::blas_int *n,
      float *                        a,
      const dealii::types::blas_int *lda,
      float *                        wr,
      float *                        wi,
      float *                        vl,
      const dealii::types::blas_int *ldvl,
      float *                        vr,
      const dealii::types::blas_int *ldvr,
      dealii::types::blas_int *      ilo,
      dealii::types::blas_int *      ihi,
      float *                        scale,
      float *                        abnrm,
      float *                        rconde,
      float *                        rcondv,
      float *                        work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  sgeevx_(balanc,
          jobvl,
          jobvr,
          sense,
          n,
          a,
          lda,
          wr,
          wi,
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
#else
  (void)balanc;
  (void)jobvl;
  (void)jobvr;
  (void)sense;
  (void)n;
  (void)a;
  (void)lda;
  (void)wr;
  (void)wi;
  (void)vl;
  (void)ldvl;
  (void)vr;
  (void)ldvr;
  (void)ilo;
  (void)ihi;
  (void)scale;
  (void)abnrm;
  (void)rconde;
  (void)rcondv;
  (void)work;
  (void)lwork;
  (void)iwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("sgeevx"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
geevx(const char *                   balanc,
      const char *                   jobvl,
      const char *                   jobvr,
      const char *                   sense,
      const dealii::types::blas_int *n,
      double *                       a,
      const dealii::types::blas_int *lda,
      double *                       wr,
      double *                       wi,
      double *                       vl,
      const dealii::types::blas_int *ldvl,
      double *                       vr,
      const dealii::types::blas_int *ldvr,
      dealii::types::blas_int *      ilo,
      dealii::types::blas_int *      ihi,
      double *                       scale,
      double *                       abnrm,
      double *                       rconde,
      double *                       rcondv,
      double *                       work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dgeevx_(balanc,
          jobvl,
          jobvr,
          sense,
          n,
          a,
          lda,
          wr,
          wi,
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
#else
  (void)balanc;
  (void)jobvl;
  (void)jobvr;
  (void)sense;
  (void)n;
  (void)a;
  (void)lda;
  (void)wr;
  (void)wi;
  (void)vl;
  (void)ldvl;
  (void)vr;
  (void)ldvr;
  (void)ilo;
  (void)ihi;
  (void)scale;
  (void)abnrm;
  (void)rconde;
  (void)rcondv;
  (void)work;
  (void)lwork;
  (void)iwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dgeevx"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
geevx(const char *                   balanc,
      const char *                   jobvl,
      const char *                   jobvr,
      const char *                   sense,
      const dealii::types::blas_int *n,
      std::complex<float> *          a,
      const dealii::types::blas_int *lda,
      std::complex<float> *          w,
      std::complex<float> *          vl,
      const dealii::types::blas_int *ldvl,
      std::complex<float> *          vr,
      const dealii::types::blas_int *ldvr,
      dealii::types::blas_int *      ilo,
      dealii::types::blas_int *      ihi,
      float *                        scale,
      float *                        abnrm,
      float *                        rconde,
      float *                        rcondv,
      std::complex<float> *          work,
      const dealii::types::blas_int *lwork,
      float *                        rwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  cgeevx_(balanc,
          jobvl,
          jobvr,
          sense,
          n,
          a,
          lda,
          w,
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
          rwork,
          info);
#else
  (void)balanc;
  (void)jobvl;
  (void)jobvr;
  (void)sense;
  (void)n;
  (void)a;
  (void)lda;
  (void)w;
  (void)vl;
  (void)ldvl;
  (void)vr;
  (void)ldvr;
  (void)ilo;
  (void)ihi;
  (void)scale;
  (void)abnrm;
  (void)rconde;
  (void)rcondv;
  (void)work;
  (void)lwork;
  (void)rwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("cgeevx"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
geevx(const char *                   balanc,
      const char *                   jobvl,
      const char *                   jobvr,
      const char *                   sense,
      const dealii::types::blas_int *n,
      std::complex<double> *         a,
      const dealii::types::blas_int *lda,
      std::complex<double> *         w,
      std::complex<double> *         vl,
      const dealii::types::blas_int *ldvl,
      std::complex<double> *         vr,
      const dealii::types::blas_int *ldvr,
      dealii::types::blas_int *      ilo,
      dealii::types::blas_int *      ihi,
      double *                       scale,
      double *                       abnrm,
      double *                       rconde,
      double *                       rcondv,
      std::complex<double> *         work,
      const dealii::types::blas_int *lwork,
      double *                       rwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  zgeevx_(balanc,
          jobvl,
          jobvr,
          sense,
          n,
          a,
          lda,
          w,
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
          rwork,
          info);
#else
  (void)balanc;
  (void)jobvl;
  (void)jobvr;
  (void)sense;
  (void)n;
  (void)a;
  (void)lda;
  (void)w;
  (void)vl;
  (void)ldvl;
  (void)vr;
  (void)ldvr;
  (void)ilo;
  (void)ihi;
  (void)scale;
  (void)abnrm;
  (void)rconde;
  (void)rcondv;
  (void)work;
  (void)lwork;
  (void)rwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("zgeevx"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5>
inline void
gelsd(const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      number1 *,
      const dealii::types::blas_int *,
      number2 *,
      const dealii::types::blas_int *,
      number3 *,
      const number4 *,
      dealii::types::blas_int *,
      number5 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
gelsd(const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      float *                        a,
      const dealii::types::blas_int *lda,
      float *                        b,
      const dealii::types::blas_int *ldb,
      float *                        s,
      const float *                  rcond,
      dealii::types::blas_int *      rank,
      float *                        work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  sgelsd_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info);
#else
  (void)m;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)s;
  (void)rcond;
  (void)rank;
  (void)work;
  (void)lwork;
  (void)iwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("sgelsd"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
gelsd(const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      double *                       a,
      const dealii::types::blas_int *lda,
      double *                       b,
      const dealii::types::blas_int *ldb,
      double *                       s,
      const double *                 rcond,
      dealii::types::blas_int *      rank,
      double *                       work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dgelsd_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info);
#else
  (void)m;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)s;
  (void)rcond;
  (void)rank;
  (void)work;
  (void)lwork;
  (void)iwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dgelsd"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
gelsd(const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      std::complex<float> *          a,
      const dealii::types::blas_int *lda,
      std::complex<float> *          b,
      const dealii::types::blas_int *ldb,
      float *                        s,
      const float *                  rcond,
      dealii::types::blas_int *      rank,
      std::complex<float> *          work,
      const dealii::types::blas_int *lwork,
      float *                        rwork,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  cgelsd_(m,
          n,
          nrhs,
          a,
          lda,
          b,
          ldb,
          s,
          rcond,
          rank,
          work,
          lwork,
          rwork,
          iwork,
          info);
#else
  (void)m;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)s;
  (void)rcond;
  (void)rank;
  (void)work;
  (void)lwork;
  (void)rwork;
  (void)iwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("cgelsd"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
gelsd(const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      std::complex<double> *         a,
      const dealii::types::blas_int *lda,
      std::complex<double> *         b,
      const dealii::types::blas_int *ldb,
      double *                       s,
      const double *                 rcond,
      dealii::types::blas_int *      rank,
      std::complex<double> *         work,
      const dealii::types::blas_int *lwork,
      double *                       rwork,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  zgelsd_(m,
          n,
          nrhs,
          a,
          lda,
          b,
          ldb,
          s,
          rcond,
          rank,
          work,
          lwork,
          rwork,
          iwork,
          info);
#else
  (void)m;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)s;
  (void)rcond;
  (void)rank;
  (void)work;
  (void)lwork;
  (void)rwork;
  (void)iwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("zgelsd"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5>
inline void
gemm(const char *,
     const char *,
     const dealii::types::blas_int *,
     const dealii::types::blas_int *,
     const dealii::types::blas_int *,
     const number1 *,
     const number2 *,
     const dealii::types::blas_int *,
     const number3 *,
     const dealii::types::blas_int *,
     const number4 *,
     number5 *,
     const dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
gemm(const char *                   transa,
     const char *                   transb,
     const dealii::types::blas_int *m,
     const dealii::types::blas_int *n,
     const dealii::types::blas_int *k,
     const float *                  alpha,
     const float *                  a,
     const dealii::types::blas_int *lda,
     const float *                  b,
     const dealii::types::blas_int *ldb,
     const float *                  beta,
     float *                        c,
     const dealii::types::blas_int *ldc)
{
#ifdef DEAL_II_WITH_LAPACK
  sgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
  (void)transa;
  (void)transb;
  (void)m;
  (void)n;
  (void)k;
  (void)alpha;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)beta;
  (void)c;
  (void)ldc;
  Assert(false, LAPACKSupport::ExcMissing("sgemm"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
gemm(const char *                   transa,
     const char *                   transb,
     const dealii::types::blas_int *m,
     const dealii::types::blas_int *n,
     const dealii::types::blas_int *k,
     const double *                 alpha,
     const double *                 a,
     const dealii::types::blas_int *lda,
     const double *                 b,
     const dealii::types::blas_int *ldb,
     const double *                 beta,
     double *                       c,
     const dealii::types::blas_int *ldc)
{
#ifdef DEAL_II_WITH_LAPACK
  dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
  (void)transa;
  (void)transb;
  (void)m;
  (void)n;
  (void)k;
  (void)alpha;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)beta;
  (void)c;
  (void)ldc;
  Assert(false, LAPACKSupport::ExcMissing("dgemm"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
gemm(const char *                   transa,
     const char *                   transb,
     const dealii::types::blas_int *m,
     const dealii::types::blas_int *n,
     const dealii::types::blas_int *k,
     const std::complex<float> *    alpha,
     const std::complex<float> *    a,
     const dealii::types::blas_int *lda,
     const std::complex<float> *    b,
     const dealii::types::blas_int *ldb,
     const std::complex<float> *    beta,
     std::complex<float> *          c,
     const dealii::types::blas_int *ldc)
{
#ifdef DEAL_II_WITH_LAPACK
  cgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
  (void)transa;
  (void)transb;
  (void)m;
  (void)n;
  (void)k;
  (void)alpha;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)beta;
  (void)c;
  (void)ldc;
  Assert(false, LAPACKSupport::ExcMissing("cgemm"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
gemm(const char *                   transa,
     const char *                   transb,
     const dealii::types::blas_int *m,
     const dealii::types::blas_int *n,
     const dealii::types::blas_int *k,
     const std::complex<double> *   alpha,
     const std::complex<double> *   a,
     const dealii::types::blas_int *lda,
     const std::complex<double> *   b,
     const dealii::types::blas_int *ldb,
     const std::complex<double> *   beta,
     std::complex<double> *         c,
     const dealii::types::blas_int *ldc)
{
#ifdef DEAL_II_WITH_LAPACK
  zgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
  (void)transa;
  (void)transb;
  (void)m;
  (void)n;
  (void)k;
  (void)alpha;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)beta;
  (void)c;
  (void)ldc;
  Assert(false, LAPACKSupport::ExcMissing("zgemm"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5>
inline void
gemv(const char *,
     const dealii::types::blas_int *,
     const dealii::types::blas_int *,
     const number1 *,
     const number2 *,
     const dealii::types::blas_int *,
     const number3 *,
     const dealii::types::blas_int *,
     const number4 *,
     number5 *,
     const dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
gemv(const char *                   trans,
     const dealii::types::blas_int *m,
     const dealii::types::blas_int *n,
     const float *                  alpha,
     const float *                  a,
     const dealii::types::blas_int *lda,
     const float *                  x,
     const dealii::types::blas_int *incx,
     const float *                  beta,
     float *                        y,
     const dealii::types::blas_int *incy)
{
#ifdef DEAL_II_WITH_LAPACK
  sgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#else
  (void)trans;
  (void)m;
  (void)n;
  (void)alpha;
  (void)a;
  (void)lda;
  (void)x;
  (void)incx;
  (void)beta;
  (void)y;
  (void)incy;
  Assert(false, LAPACKSupport::ExcMissing("sgemv"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
gemv(const char *                   trans,
     const dealii::types::blas_int *m,
     const dealii::types::blas_int *n,
     const double *                 alpha,
     const double *                 a,
     const dealii::types::blas_int *lda,
     const double *                 x,
     const dealii::types::blas_int *incx,
     const double *                 beta,
     double *                       y,
     const dealii::types::blas_int *incy)
{
#ifdef DEAL_II_WITH_LAPACK
  dgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#else
  (void)trans;
  (void)m;
  (void)n;
  (void)alpha;
  (void)a;
  (void)lda;
  (void)x;
  (void)incx;
  (void)beta;
  (void)y;
  (void)incy;
  Assert(false, LAPACKSupport::ExcMissing("dgemv"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
gemv(const char *                   trans,
     const dealii::types::blas_int *m,
     const dealii::types::blas_int *n,
     const std::complex<float> *    alpha,
     const std::complex<float> *    a,
     const dealii::types::blas_int *lda,
     const std::complex<float> *    x,
     const dealii::types::blas_int *incx,
     const std::complex<float> *    beta,
     std::complex<float> *          y,
     const dealii::types::blas_int *incy)
{
#ifdef DEAL_II_WITH_LAPACK
  cgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#else
  (void)trans;
  (void)m;
  (void)n;
  (void)alpha;
  (void)a;
  (void)lda;
  (void)x;
  (void)incx;
  (void)beta;
  (void)y;
  (void)incy;
  Assert(false, LAPACKSupport::ExcMissing("cgemv"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
gemv(const char *                   trans,
     const dealii::types::blas_int *m,
     const dealii::types::blas_int *n,
     const std::complex<double> *   alpha,
     const std::complex<double> *   a,
     const dealii::types::blas_int *lda,
     const std::complex<double> *   x,
     const dealii::types::blas_int *incx,
     const std::complex<double> *   beta,
     std::complex<double> *         y,
     const dealii::types::blas_int *incy)
{
#ifdef DEAL_II_WITH_LAPACK
  zgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#else
  (void)trans;
  (void)m;
  (void)n;
  (void)alpha;
  (void)a;
  (void)lda;
  (void)x;
  (void)incx;
  (void)beta;
  (void)y;
  (void)incy;
  Assert(false, LAPACKSupport::ExcMissing("zgemv"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1, typename number2, typename number3>
inline void
geqrf(const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      number1 *,
      const dealii::types::blas_int *,
      number2 *,
      number3 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
geqrf(const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      float *                        a,
      const dealii::types::blas_int *lda,
      float *                        tau,
      float *                        work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  sgeqrf_(m, n, a, lda, tau, work, lwork, info);
#else
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)tau;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("sgeqrf"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
geqrf(const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      double *                       a,
      const dealii::types::blas_int *lda,
      double *                       tau,
      double *                       work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dgeqrf_(m, n, a, lda, tau, work, lwork, info);
#else
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)tau;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dgeqrf"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
geqrf(const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      std::complex<float> *          a,
      const dealii::types::blas_int *lda,
      std::complex<float> *          tau,
      std::complex<float> *          work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  cgeqrf_(m, n, a, lda, tau, work, lwork, info);
#else
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)tau;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("cgeqrf"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
geqrf(const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      std::complex<double> *         a,
      const dealii::types::blas_int *lda,
      std::complex<double> *         tau,
      std::complex<double> *         work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  zgeqrf_(m, n, a, lda, tau, work, lwork, info);
#else
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)tau;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("zgeqrf"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5>
inline void
gesdd(const char *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      number1 *,
      const dealii::types::blas_int *,
      number2 *,
      number3 *,
      const dealii::types::blas_int *,
      number4 *,
      const dealii::types::blas_int *,
      number5 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
gesdd(const char *                   jobz,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      float *                        a,
      const dealii::types::blas_int *lda,
      float *                        s,
      float *                        u,
      const dealii::types::blas_int *ldu,
      float *                        vt,
      const dealii::types::blas_int *ldvt,
      float *                        work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  sgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
#else
  (void)jobz;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)s;
  (void)u;
  (void)ldu;
  (void)vt;
  (void)ldvt;
  (void)work;
  (void)lwork;
  (void)iwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("sgesdd"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
gesdd(const char *                   jobz,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      double *                       a,
      const dealii::types::blas_int *lda,
      double *                       s,
      double *                       u,
      const dealii::types::blas_int *ldu,
      double *                       vt,
      const dealii::types::blas_int *ldvt,
      double *                       work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
#else
  (void)jobz;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)s;
  (void)u;
  (void)ldu;
  (void)vt;
  (void)ldvt;
  (void)work;
  (void)lwork;
  (void)iwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dgesdd"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
gesdd(const char *                   jobz,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      std::complex<float> *          a,
      const dealii::types::blas_int *lda,
      float *                        s,
      std::complex<float> *          u,
      const dealii::types::blas_int *ldu,
      std::complex<float> *          vt,
      const dealii::types::blas_int *ldvt,
      std::complex<float> *          work,
      const dealii::types::blas_int *lwork,
      float *                        rwork,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  cgesdd_(
    jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
#else
  (void)jobz;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)s;
  (void)u;
  (void)ldu;
  (void)vt;
  (void)ldvt;
  (void)work;
  (void)lwork;
  (void)rwork;
  (void)iwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("cgesdd"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
gesdd(const char *                   jobz,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      std::complex<double> *         a,
      const dealii::types::blas_int *lda,
      double *                       s,
      std::complex<double> *         u,
      const dealii::types::blas_int *ldu,
      std::complex<double> *         vt,
      const dealii::types::blas_int *ldvt,
      std::complex<double> *         work,
      const dealii::types::blas_int *lwork,
      double *                       rwork,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  zgesdd_(
    jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
#else
  (void)jobz;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)s;
  (void)u;
  (void)ldu;
  (void)vt;
  (void)ldvt;
  (void)work;
  (void)lwork;
  (void)rwork;
  (void)iwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("zgesdd"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5>
inline void
gesvd(const char *,
      const char *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      number1 *,
      const dealii::types::blas_int *,
      number2 *,
      number3 *,
      const dealii::types::blas_int *,
      number4 *,
      const dealii::types::blas_int *,
      number5 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
gesvd(const char *                   jobu,
      const char *                   jobvt,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      float *                        a,
      const dealii::types::blas_int *lda,
      float *                        s,
      float *                        u,
      const dealii::types::blas_int *ldu,
      float *                        vt,
      const dealii::types::blas_int *ldvt,
      float *                        work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  sgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
#else
  (void)jobu;
  (void)jobvt;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)s;
  (void)u;
  (void)ldu;
  (void)vt;
  (void)ldvt;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("sgesvd"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
gesvd(const char *                   jobu,
      const char *                   jobvt,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      double *                       a,
      const dealii::types::blas_int *lda,
      double *                       s,
      double *                       u,
      const dealii::types::blas_int *ldu,
      double *                       vt,
      const dealii::types::blas_int *ldvt,
      double *                       work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
#else
  (void)jobu;
  (void)jobvt;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)s;
  (void)u;
  (void)ldu;
  (void)vt;
  (void)ldvt;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dgesvd"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
gesvd(const char *                   jobu,
      const char *                   jobvt,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      std::complex<float> *          a,
      const dealii::types::blas_int *lda,
      float *                        s,
      std::complex<float> *          u,
      const dealii::types::blas_int *ldu,
      std::complex<float> *          vt,
      const dealii::types::blas_int *ldvt,
      std::complex<float> *          work,
      const dealii::types::blas_int *lwork,
      float *                        rwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  cgesvd_(
    jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
#else
  (void)jobu;
  (void)jobvt;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)s;
  (void)u;
  (void)ldu;
  (void)vt;
  (void)ldvt;
  (void)work;
  (void)lwork;
  (void)rwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("cgesvd"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
gesvd(const char *                   jobu,
      const char *                   jobvt,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      std::complex<double> *         a,
      const dealii::types::blas_int *lda,
      double *                       s,
      std::complex<double> *         u,
      const dealii::types::blas_int *ldu,
      std::complex<double> *         vt,
      const dealii::types::blas_int *ldvt,
      std::complex<double> *         work,
      const dealii::types::blas_int *lwork,
      double *                       rwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  zgesvd_(
    jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
#else
  (void)jobu;
  (void)jobvt;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)s;
  (void)u;
  (void)ldu;
  (void)vt;
  (void)ldvt;
  (void)work;
  (void)lwork;
  (void)rwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("zgesvd"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1>
inline void
getrf(const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      number1 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
getrf(const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      float *                        a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      ipiv,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  sgetrf_(m, n, a, lda, ipiv, info);
#else
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)ipiv;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("sgetrf"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
getrf(const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      double *                       a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      ipiv,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dgetrf_(m, n, a, lda, ipiv, info);
#else
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)ipiv;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dgetrf"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
getrf(const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      std::complex<float> *          a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      ipiv,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  cgetrf_(m, n, a, lda, ipiv, info);
#else
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)ipiv;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("cgetrf"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
getrf(const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      std::complex<double> *         a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      ipiv,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  zgetrf_(m, n, a, lda, ipiv, info);
#else
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)ipiv;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("zgetrf"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1, typename number2>
inline void
getri(const dealii::types::blas_int *,
      number1 *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      number2 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
getri(const dealii::types::blas_int *n,
      float *                        a,
      const dealii::types::blas_int *lda,
      const dealii::types::blas_int *ipiv,
      float *                        work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  sgetri_(n, a, lda, ipiv, work, lwork, info);
#else
  (void)n;
  (void)a;
  (void)lda;
  (void)ipiv;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("sgetri"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
getri(const dealii::types::blas_int *n,
      double *                       a,
      const dealii::types::blas_int *lda,
      const dealii::types::blas_int *ipiv,
      double *                       work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dgetri_(n, a, lda, ipiv, work, lwork, info);
#else
  (void)n;
  (void)a;
  (void)lda;
  (void)ipiv;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dgetri"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
getri(const dealii::types::blas_int *n,
      std::complex<float> *          a,
      const dealii::types::blas_int *lda,
      const dealii::types::blas_int *ipiv,
      std::complex<float> *          work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  cgetri_(n, a, lda, ipiv, work, lwork, info);
#else
  (void)n;
  (void)a;
  (void)lda;
  (void)ipiv;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("cgetri"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
getri(const dealii::types::blas_int *n,
      std::complex<double> *         a,
      const dealii::types::blas_int *lda,
      const dealii::types::blas_int *ipiv,
      std::complex<double> *         work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  zgetri_(n, a, lda, ipiv, work, lwork, info);
#else
  (void)n;
  (void)a;
  (void)lda;
  (void)ipiv;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("zgetri"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1, typename number2>
inline void
getrs(const char *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      const number1 *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      number2 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
getrs(const char *                   trans,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      const float *                  a,
      const dealii::types::blas_int *lda,
      const dealii::types::blas_int *ipiv,
      float *                        b,
      const dealii::types::blas_int *ldb,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  sgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
#else
  (void)trans;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)ipiv;
  (void)b;
  (void)ldb;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("sgetrs"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
getrs(const char *                   trans,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      const double *                 a,
      const dealii::types::blas_int *lda,
      const dealii::types::blas_int *ipiv,
      double *                       b,
      const dealii::types::blas_int *ldb,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
#else
  (void)trans;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)ipiv;
  (void)b;
  (void)ldb;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dgetrs"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
getrs(const char *                   trans,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      const std::complex<float> *    a,
      const dealii::types::blas_int *lda,
      const dealii::types::blas_int *ipiv,
      std::complex<float> *          b,
      const dealii::types::blas_int *ldb,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  cgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
#else
  (void)trans;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)ipiv;
  (void)b;
  (void)ldb;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("cgetrs"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
getrs(const char *                   trans,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      const std::complex<double> *   a,
      const dealii::types::blas_int *lda,
      const dealii::types::blas_int *ipiv,
      std::complex<double> *         b,
      const dealii::types::blas_int *ldb,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  zgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
#else
  (void)trans;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)ipiv;
  (void)b;
  (void)ldb;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("zgetrs"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1>
inline number1
lamch(const char *)
{
  Assert(false, ExcNotImplemented());
  return number1();
}



template <>
inline float
lamch(const char *cmach)
{
#ifdef DEAL_II_WITH_LAPACK
  return slamch_(cmach);
#else
  (void)cmach;
  Assert(false, LAPACKSupport::ExcMissing("slamch"));
  return float();
#endif // DEAL_II_WITH_LAPACK
}



template <>
inline double
lamch(const char *cmach)
{
#ifdef DEAL_II_WITH_LAPACK
  return dlamch_(cmach);
#else
  (void)cmach;
  Assert(false, LAPACKSupport::ExcMissing("dlamch"));
  return double();
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1, typename number2>
inline number1
lange(const char *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      const number1 *,
      const dealii::types::blas_int *,
      number2 *)
{
  Assert(false, ExcNotImplemented());
  return number1();
}



inline float
lange(const char *                   norm,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      const float *                  a,
      const dealii::types::blas_int *lda,
      float *                        work)
{
#ifdef DEAL_II_WITH_LAPACK
  return slange_(norm, m, n, a, lda, work);
#else
  (void)norm;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)work;
  Assert(false, LAPACKSupport::ExcMissing("slange"));
  return float();
#endif // DEAL_II_WITH_LAPACK
}



inline double
lange(const char *                   norm,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      const double *                 a,
      const dealii::types::blas_int *lda,
      double *                       work)
{
#ifdef DEAL_II_WITH_LAPACK
  return dlange_(norm, m, n, a, lda, work);
#else
  (void)norm;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)work;
  Assert(false, LAPACKSupport::ExcMissing("dlange"));
  return double();
#endif // DEAL_II_WITH_LAPACK
}



inline float
lange(const char *                   norm,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      const std::complex<float> *    a,
      const dealii::types::blas_int *lda,
      float *                        work)
{
#ifdef DEAL_II_WITH_LAPACK
  return clange_(norm, m, n, a, lda, work);
#else
  (void)norm;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)work;
  Assert(false, LAPACKSupport::ExcMissing("clange"));
  return float();
#endif // DEAL_II_WITH_LAPACK
}



inline double
lange(const char *                   norm,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      const std::complex<double> *   a,
      const dealii::types::blas_int *lda,
      double *                       work)
{
#ifdef DEAL_II_WITH_LAPACK
  return zlange_(norm, m, n, a, lda, work);
#else
  (void)norm;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)work;
  Assert(false, LAPACKSupport::ExcMissing("zlange"));
  return double();
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1, typename number2>
inline number1
lansy(const char *,
      const char *,
      const dealii::types::blas_int *,
      const number1 *,
      const dealii::types::blas_int *,
      number2 *)
{
  Assert(false, ExcNotImplemented());
  return number1();
}



inline float
lansy(const char *                   norm,
      const char *                   uplo,
      const dealii::types::blas_int *n,
      const float *                  a,
      const dealii::types::blas_int *lda,
      float *                        work)
{
#ifdef DEAL_II_WITH_LAPACK
  return slansy_(norm, uplo, n, a, lda, work);
#else
  (void)norm;
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)work;
  Assert(false, LAPACKSupport::ExcMissing("slansy"));
  return float();
#endif // DEAL_II_WITH_LAPACK
}



inline double
lansy(const char *                   norm,
      const char *                   uplo,
      const dealii::types::blas_int *n,
      const double *                 a,
      const dealii::types::blas_int *lda,
      double *                       work)
{
#ifdef DEAL_II_WITH_LAPACK
  return dlansy_(norm, uplo, n, a, lda, work);
#else
  (void)norm;
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)work;
  Assert(false, LAPACKSupport::ExcMissing("dlansy"));
  return double();
#endif // DEAL_II_WITH_LAPACK
}



inline float
lansy(const char *                   norm,
      const char *                   uplo,
      const dealii::types::blas_int *n,
      const std::complex<float> *    a,
      const dealii::types::blas_int *lda,
      float *                        work)
{
#ifdef DEAL_II_WITH_LAPACK
  return clansy_(norm, uplo, n, a, lda, work);
#else
  (void)norm;
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)work;
  Assert(false, LAPACKSupport::ExcMissing("clansy"));
  return float();
#endif // DEAL_II_WITH_LAPACK
}



inline double
lansy(const char *                   norm,
      const char *                   uplo,
      const dealii::types::blas_int *n,
      const std::complex<double> *   a,
      const dealii::types::blas_int *lda,
      double *                       work)
{
#ifdef DEAL_II_WITH_LAPACK
  return zlansy_(norm, uplo, n, a, lda, work);
#else
  (void)norm;
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)work;
  Assert(false, LAPACKSupport::ExcMissing("zlansy"));
  return double();
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1, typename number2, typename number3>
inline void
lascl(const char *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      const number1 *,
      const number2 *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      number3 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
lascl(const char *                   type,
      const dealii::types::blas_int *kl,
      const dealii::types::blas_int *ku,
      const float *                  cfrom,
      const float *                  cto,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      float *                        a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  slascl_(type, kl, ku, cfrom, cto, m, n, a, lda, info);
#else
  (void)type;
  (void)kl;
  (void)ku;
  (void)cfrom;
  (void)cto;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("slascl"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
lascl(const char *                   type,
      const dealii::types::blas_int *kl,
      const dealii::types::blas_int *ku,
      const double *                 cfrom,
      const double *                 cto,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      double *                       a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dlascl_(type, kl, ku, cfrom, cto, m, n, a, lda, info);
#else
  (void)type;
  (void)kl;
  (void)ku;
  (void)cfrom;
  (void)cto;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dlascl"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
lascl(const char *                   type,
      const dealii::types::blas_int *kl,
      const dealii::types::blas_int *ku,
      const float *                  cfrom,
      const float *                  cto,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      std::complex<float> *          a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  clascl_(type, kl, ku, cfrom, cto, m, n, a, lda, info);
#else
  (void)type;
  (void)kl;
  (void)ku;
  (void)cfrom;
  (void)cto;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("clascl"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
lascl(const char *                   type,
      const dealii::types::blas_int *kl,
      const dealii::types::blas_int *ku,
      const double *                 cfrom,
      const double *                 cto,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      std::complex<double> *         a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  zlascl_(type, kl, ku, cfrom, cto, m, n, a, lda, info);
#else
  (void)type;
  (void)kl;
  (void)ku;
  (void)cfrom;
  (void)cto;
  (void)m;
  (void)n;
  (void)a;
  (void)lda;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("zlascl"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1, typename number2, typename number3>
inline void
orgqr(const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      number1 *,
      const dealii::types::blas_int *,
      const number2 *,
      number3 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
orgqr(const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *k,
      float *                        a,
      const dealii::types::blas_int *lda,
      const float *                  tau,
      float *                        work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  sorgqr_(m, n, k, a, lda, tau, work, lwork, info);
#else
  (void)m;
  (void)n;
  (void)k;
  (void)a;
  (void)lda;
  (void)tau;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("sorgqr"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
orgqr(const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *k,
      double *                       a,
      const dealii::types::blas_int *lda,
      const double *                 tau,
      double *                       work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dorgqr_(m, n, k, a, lda, tau, work, lwork, info);
#else
  (void)m;
  (void)n;
  (void)k;
  (void)a;
  (void)lda;
  (void)tau;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dorgqr"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1,
          typename number2,
          typename number3,
          typename number4>
inline void
ormqr(const char *,
      const char *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      const number1 *,
      const dealii::types::blas_int *,
      const number2 *,
      number3 *,
      const dealii::types::blas_int *,
      number4 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
ormqr(const char *                   side,
      const char *                   trans,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *k,
      const float *                  a,
      const dealii::types::blas_int *lda,
      const float *                  tau,
      float *                        c,
      const dealii::types::blas_int *ldc,
      float *                        work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  sormqr_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
#else
  (void)side;
  (void)trans;
  (void)m;
  (void)n;
  (void)k;
  (void)a;
  (void)lda;
  (void)tau;
  (void)c;
  (void)ldc;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("sormqr"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
ormqr(const char *                   side,
      const char *                   trans,
      const dealii::types::blas_int *m,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *k,
      const double *                 a,
      const dealii::types::blas_int *lda,
      const double *                 tau,
      double *                       c,
      const dealii::types::blas_int *ldc,
      double *                       work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dormqr_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
#else
  (void)side;
  (void)trans;
  (void)m;
  (void)n;
  (void)k;
  (void)a;
  (void)lda;
  (void)tau;
  (void)c;
  (void)ldc;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dormqr"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1,
          typename number2,
          typename number3,
          typename number4>
inline void
pocon(const char *,
      const dealii::types::blas_int *,
      const number1 *,
      const dealii::types::blas_int *,
      const number2 *,
      number3 *,
      number4 *,
      dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
pocon(const char *                   uplo,
      const dealii::types::blas_int *n,
      const float *                  a,
      const dealii::types::blas_int *lda,
      const float *                  anorm,
      float *                        rcond,
      float *                        work,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  spocon_(uplo, n, a, lda, anorm, rcond, work, iwork, info);
#else
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)anorm;
  (void)rcond;
  (void)work;
  (void)iwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("spocon"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
pocon(const char *                   uplo,
      const dealii::types::blas_int *n,
      const double *                 a,
      const dealii::types::blas_int *lda,
      const double *                 anorm,
      double *                       rcond,
      double *                       work,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dpocon_(uplo, n, a, lda, anorm, rcond, work, iwork, info);
#else
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)anorm;
  (void)rcond;
  (void)work;
  (void)iwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dpocon"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
pocon(const char *                   uplo,
      const dealii::types::blas_int *n,
      const std::complex<float> *    a,
      const dealii::types::blas_int *lda,
      const float *                  anorm,
      float *                        rcond,
      std::complex<float> *          work,
      float *                        rwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  cpocon_(uplo, n, a, lda, anorm, rcond, work, rwork, info);
#else
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)anorm;
  (void)rcond;
  (void)work;
  (void)rwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("cpocon"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
pocon(const char *                   uplo,
      const dealii::types::blas_int *n,
      const std::complex<double> *   a,
      const dealii::types::blas_int *lda,
      const double *                 anorm,
      double *                       rcond,
      std::complex<double> *         work,
      double *                       rwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  zpocon_(uplo, n, a, lda, anorm, rcond, work, rwork, info);
#else
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)anorm;
  (void)rcond;
  (void)work;
  (void)rwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("zpocon"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1>
inline void
potrf(const char *,
      const dealii::types::blas_int *,
      number1 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
potrf(const char *                   uplo,
      const dealii::types::blas_int *n,
      float *                        a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  spotrf_(uplo, n, a, lda, info);
#else
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("spotrf"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
potrf(const char *                   uplo,
      const dealii::types::blas_int *n,
      double *                       a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dpotrf_(uplo, n, a, lda, info);
#else
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dpotrf"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
potrf(const char *                   uplo,
      const dealii::types::blas_int *n,
      std::complex<float> *          a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  cpotrf_(uplo, n, a, lda, info);
#else
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("cpotrf"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
potrf(const char *                   uplo,
      const dealii::types::blas_int *n,
      std::complex<double> *         a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  zpotrf_(uplo, n, a, lda, info);
#else
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("zpotrf"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1>
inline void
potri(const char *,
      const dealii::types::blas_int *,
      number1 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
potri(const char *                   uplo,
      const dealii::types::blas_int *n,
      float *                        a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  spotri_(uplo, n, a, lda, info);
#else
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("spotri"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
potri(const char *                   uplo,
      const dealii::types::blas_int *n,
      double *                       a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dpotri_(uplo, n, a, lda, info);
#else
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dpotri"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
potri(const char *                   uplo,
      const dealii::types::blas_int *n,
      std::complex<float> *          a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  cpotri_(uplo, n, a, lda, info);
#else
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("cpotri"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
potri(const char *                   uplo,
      const dealii::types::blas_int *n,
      std::complex<double> *         a,
      const dealii::types::blas_int *lda,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  zpotri_(uplo, n, a, lda, info);
#else
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("zpotri"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1, typename number2>
inline void
potrs(const char *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      const number1 *,
      const dealii::types::blas_int *,
      number2 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
potrs(const char *                   uplo,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      const float *                  a,
      const dealii::types::blas_int *lda,
      float *                        b,
      const dealii::types::blas_int *ldb,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  spotrs_(uplo, n, nrhs, a, lda, b, ldb, info);
#else
  (void)uplo;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("spotrs"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
potrs(const char *                   uplo,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      const double *                 a,
      const dealii::types::blas_int *lda,
      double *                       b,
      const dealii::types::blas_int *ldb,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dpotrs_(uplo, n, nrhs, a, lda, b, ldb, info);
#else
  (void)uplo;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dpotrs"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
potrs(const char *                   uplo,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      const std::complex<float> *    a,
      const dealii::types::blas_int *lda,
      std::complex<float> *          b,
      const dealii::types::blas_int *ldb,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  cpotrs_(uplo, n, nrhs, a, lda, b, ldb, info);
#else
  (void)uplo;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("cpotrs"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
potrs(const char *                   uplo,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      const std::complex<double> *   a,
      const dealii::types::blas_int *lda,
      std::complex<double> *         b,
      const dealii::types::blas_int *ldb,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  zpotrs_(uplo, n, nrhs, a, lda, b, ldb, info);
#else
  (void)uplo;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("zpotrs"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1,
          typename number2,
          typename number3,
          typename number4>
inline void
stev(const char *,
     const dealii::types::blas_int *,
     number1 *,
     number2 *,
     number3 *,
     const dealii::types::blas_int *,
     number4 *,
     dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
stev(const char *                   jobz,
     const dealii::types::blas_int *n,
     float *                        d,
     float *                        e,
     float *                        z,
     const dealii::types::blas_int *ldz,
     float *                        work,
     dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  sstev_(jobz, n, d, e, z, ldz, work, info);
#else
  (void)jobz;
  (void)n;
  (void)d;
  (void)e;
  (void)z;
  (void)ldz;
  (void)work;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("sstev"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
stev(const char *                   jobz,
     const dealii::types::blas_int *n,
     double *                       d,
     double *                       e,
     double *                       z,
     const dealii::types::blas_int *ldz,
     double *                       work,
     dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dstev_(jobz, n, d, e, z, ldz, work, info);
#else
  (void)jobz;
  (void)n;
  (void)d;
  (void)e;
  (void)z;
  (void)ldz;
  (void)work;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dstev"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1, typename number2, typename number3>
inline void
syev(const char *,
     const char *,
     const dealii::types::blas_int *,
     number1 *,
     const dealii::types::blas_int *,
     number2 *,
     number3 *,
     const dealii::types::blas_int *,
     dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
syev(const char *                   jobz,
     const char *                   uplo,
     const dealii::types::blas_int *n,
     float *                        a,
     const dealii::types::blas_int *lda,
     float *                        w,
     float *                        work,
     const dealii::types::blas_int *lwork,
     dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  ssyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
#else
  (void)jobz;
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)w;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("ssyev"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
syev(const char *                   jobz,
     const char *                   uplo,
     const dealii::types::blas_int *n,
     double *                       a,
     const dealii::types::blas_int *lda,
     double *                       w,
     double *                       work,
     const dealii::types::blas_int *lwork,
     dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
#else
  (void)jobz;
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)w;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dsyev"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5,
          typename number6,
          typename number7>
inline void
syevr(const char *,
      const char *,
      const char *,
      const dealii::types::blas_int *,
      number1 *,
      const dealii::types::blas_int *,
      const number2 *,
      const number3 *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      const number4 *,
      dealii::types::blas_int *,
      number5 *,
      number6 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *,
      number7 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
syevr(const char *                   jobz,
      const char *                   range,
      const char *                   uplo,
      const dealii::types::blas_int *n,
      float *                        a,
      const dealii::types::blas_int *lda,
      const float *                  vl,
      const float *                  vu,
      const dealii::types::blas_int *il,
      const dealii::types::blas_int *iu,
      const float *                  abstol,
      dealii::types::blas_int *      m,
      float *                        w,
      float *                        z,
      const dealii::types::blas_int *ldz,
      dealii::types::blas_int *      isuppz,
      float *                        work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      iwork,
      const dealii::types::blas_int *liwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  // Netlib and Atlas Lapack perform floating point tests (e.g. divide-by-zero)
  // within calls to some functions, which cause floating point exceptions to be
  // thrown (at least in debug mode). Therefore, we wrap the calls into the
  // following code to suppress the exception.
#  ifdef DEAL_II_HAVE_FP_EXCEPTIONS
  fenv_t fp_exceptions;
  feholdexcept(&fp_exceptions);
#  endif // DEAL_II_HAVE_FP_EXCEPTIONS
  ssyevr_(jobz,
          range,
          uplo,
          n,
          a,
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
#  endif // DEAL_II_HAVE_FP_EXCEPTIONS
#else
  (void)jobz;
  (void)range;
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)vl;
  (void)vu;
  (void)il;
  (void)iu;
  (void)abstol;
  (void)m;
  (void)w;
  (void)z;
  (void)ldz;
  (void)isuppz;
  (void)work;
  (void)lwork;
  (void)iwork;
  (void)liwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("ssyevr"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
syevr(const char *                   jobz,
      const char *                   range,
      const char *                   uplo,
      const dealii::types::blas_int *n,
      double *                       a,
      const dealii::types::blas_int *lda,
      const double *                 vl,
      const double *                 vu,
      const dealii::types::blas_int *il,
      const dealii::types::blas_int *iu,
      const double *                 abstol,
      dealii::types::blas_int *      m,
      double *                       w,
      double *                       z,
      const dealii::types::blas_int *ldz,
      dealii::types::blas_int *      isuppz,
      double *                       work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      iwork,
      const dealii::types::blas_int *liwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  // Netlib and Atlas Lapack perform floating point tests (e.g. divide-by-zero)
  // within calls to some functions, which cause floating point exceptions to be
  // thrown (at least in debug mode). Therefore, we wrap the calls into the
  // following code to suppress the exception.
#  ifdef DEAL_II_HAVE_FP_EXCEPTIONS
  fenv_t fp_exceptions;
  feholdexcept(&fp_exceptions);
#  endif // DEAL_II_HAVE_FP_EXCEPTIONS
  dsyevr_(jobz,
          range,
          uplo,
          n,
          a,
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
#  endif // DEAL_II_HAVE_FP_EXCEPTIONS
#else
  (void)jobz;
  (void)range;
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)vl;
  (void)vu;
  (void)il;
  (void)iu;
  (void)abstol;
  (void)m;
  (void)w;
  (void)z;
  (void)ldz;
  (void)isuppz;
  (void)work;
  (void)lwork;
  (void)iwork;
  (void)liwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dsyevr"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5,
          typename number6,
          typename number7>
inline void
syevx(const char *,
      const char *,
      const char *,
      const dealii::types::blas_int *,
      number1 *,
      const dealii::types::blas_int *,
      const number2 *,
      const number3 *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      const number4 *,
      dealii::types::blas_int *,
      number5 *,
      number6 *,
      const dealii::types::blas_int *,
      number7 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *,
      dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}


inline void
syevx(const char *                   jobz,
      const char *                   range,
      const char *                   uplo,
      const dealii::types::blas_int *n,
      float *                        a,
      const dealii::types::blas_int *lda,
      const float *                  vl,
      const float *                  vu,
      const dealii::types::blas_int *il,
      const dealii::types::blas_int *iu,
      const float *                  abstol,
      dealii::types::blas_int *      m,
      float *                        w,
      float *                        z,
      const dealii::types::blas_int *ldz,
      float *                        work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      ifail,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  ssyevx_(jobz,
          range,
          uplo,
          n,
          a,
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
#else
  (void)jobz;
  (void)range;
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)vl;
  (void)vu;
  (void)il;
  (void)iu;
  (void)abstol;
  (void)m;
  (void)w;
  (void)z;
  (void)ldz;
  (void)work;
  (void)lwork;
  (void)iwork;
  (void)ifail;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("ssyevx"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
syevx(const char *                   jobz,
      const char *                   range,
      const char *                   uplo,
      const dealii::types::blas_int *n,
      double *                       a,
      const dealii::types::blas_int *lda,
      const double *                 vl,
      const double *                 vu,
      const dealii::types::blas_int *il,
      const dealii::types::blas_int *iu,
      const double *                 abstol,
      dealii::types::blas_int *      m,
      double *                       w,
      double *                       z,
      const dealii::types::blas_int *ldz,
      double *                       work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      ifail,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dsyevx_(jobz,
          range,
          uplo,
          n,
          a,
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
#else
  (void)jobz;
  (void)range;
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)vl;
  (void)vu;
  (void)il;
  (void)iu;
  (void)abstol;
  (void)m;
  (void)w;
  (void)z;
  (void)ldz;
  (void)work;
  (void)lwork;
  (void)iwork;
  (void)ifail;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dsyevx"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1,
          typename number2,
          typename number3,
          typename number4>
inline void
sygv(const dealii::types::blas_int *,
     const char *,
     const char *,
     const dealii::types::blas_int *,
     number1 *,
     const dealii::types::blas_int *,
     number2 *,
     const dealii::types::blas_int *,
     number3 *,
     number4 *,
     const dealii::types::blas_int *,
     dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
sygv(const dealii::types::blas_int *itype,
     const char *                   jobz,
     const char *                   uplo,
     const dealii::types::blas_int *n,
     float *                        a,
     const dealii::types::blas_int *lda,
     float *                        b,
     const dealii::types::blas_int *ldb,
     float *                        w,
     float *                        work,
     const dealii::types::blas_int *lwork,
     dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  ssygv_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
#else
  (void)itype;
  (void)jobz;
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)w;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("ssygv"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
sygv(const dealii::types::blas_int *itype,
     const char *                   jobz,
     const char *                   uplo,
     const dealii::types::blas_int *n,
     double *                       a,
     const dealii::types::blas_int *lda,
     double *                       b,
     const dealii::types::blas_int *ldb,
     double *                       w,
     double *                       work,
     const dealii::types::blas_int *lwork,
     dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dsygv_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
#else
  (void)itype;
  (void)jobz;
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)w;
  (void)work;
  (void)lwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dsygv"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1,
          typename number2,
          typename number3,
          typename number4,
          typename number5,
          typename number6,
          typename number7,
          typename number8>
inline void
sygvx(const dealii::types::blas_int *,
      const char *,
      const char *,
      const char *,
      const dealii::types::blas_int *,
      number1 *,
      const dealii::types::blas_int *,
      number2 *,
      const dealii::types::blas_int *,
      const number3 *,
      const number4 *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      const number5 *,
      dealii::types::blas_int *,
      number6 *,
      number7 *,
      const dealii::types::blas_int *,
      number8 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *,
      dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
sygvx(const dealii::types::blas_int *itype,
      const char *                   jobz,
      const char *                   range,
      const char *                   uplo,
      const dealii::types::blas_int *n,
      float *                        a,
      const dealii::types::blas_int *lda,
      float *                        b,
      const dealii::types::blas_int *ldb,
      const float *                  vl,
      const float *                  vu,
      const dealii::types::blas_int *il,
      const dealii::types::blas_int *iu,
      const float *                  abstol,
      dealii::types::blas_int *      m,
      float *                        w,
      float *                        z,
      const dealii::types::blas_int *ldz,
      float *                        work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      ifail,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  ssygvx_(itype,
          jobz,
          range,
          uplo,
          n,
          a,
          lda,
          b,
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
#else
  (void)itype;
  (void)jobz;
  (void)range;
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)vl;
  (void)vu;
  (void)il;
  (void)iu;
  (void)abstol;
  (void)m;
  (void)w;
  (void)z;
  (void)ldz;
  (void)work;
  (void)lwork;
  (void)iwork;
  (void)ifail;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("ssygvx"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
sygvx(const dealii::types::blas_int *itype,
      const char *                   jobz,
      const char *                   range,
      const char *                   uplo,
      const dealii::types::blas_int *n,
      double *                       a,
      const dealii::types::blas_int *lda,
      double *                       b,
      const dealii::types::blas_int *ldb,
      const double *                 vl,
      const double *                 vu,
      const dealii::types::blas_int *il,
      const dealii::types::blas_int *iu,
      const double *                 abstol,
      dealii::types::blas_int *      m,
      double *                       w,
      double *                       z,
      const dealii::types::blas_int *ldz,
      double *                       work,
      const dealii::types::blas_int *lwork,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      ifail,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dsygvx_(itype,
          jobz,
          range,
          uplo,
          n,
          a,
          lda,
          b,
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
#else
  (void)itype;
  (void)jobz;
  (void)range;
  (void)uplo;
  (void)n;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)vl;
  (void)vu;
  (void)il;
  (void)iu;
  (void)abstol;
  (void)m;
  (void)w;
  (void)z;
  (void)ldz;
  (void)work;
  (void)lwork;
  (void)iwork;
  (void)ifail;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dsygvx"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1, typename number2, typename number3>
inline void
syr(const char *,
    const dealii::types::blas_int *,
    const number1 *,
    const number2 *,
    const dealii::types::blas_int *,
    number3 *,
    const dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
syr(const char *                   uplo,
    const dealii::types::blas_int *n,
    const float *                  alpha,
    const float *                  x,
    const dealii::types::blas_int *incx,
    float *                        a,
    const dealii::types::blas_int *lda)
{
#ifdef DEAL_II_WITH_LAPACK
  ssyr_(uplo, n, alpha, x, incx, a, lda);
#else
  (void)uplo;
  (void)n;
  (void)alpha;
  (void)x;
  (void)incx;
  (void)a;
  (void)lda;
  Assert(false, LAPACKSupport::ExcMissing("ssyr"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
syr(const char *                   uplo,
    const dealii::types::blas_int *n,
    const double *                 alpha,
    const double *                 x,
    const dealii::types::blas_int *incx,
    double *                       a,
    const dealii::types::blas_int *lda)
{
#ifdef DEAL_II_WITH_LAPACK
  dsyr_(uplo, n, alpha, x, incx, a, lda);
#else
  (void)uplo;
  (void)n;
  (void)alpha;
  (void)x;
  (void)incx;
  (void)a;
  (void)lda;
  Assert(false, LAPACKSupport::ExcMissing("dsyr"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1,
          typename number2,
          typename number3,
          typename number4>
inline void
syrk(const char *,
     const char *,
     const dealii::types::blas_int *,
     const dealii::types::blas_int *,
     const number1 *,
     const number2 *,
     const dealii::types::blas_int *,
     const number3 *,
     number4 *,
     const dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
syrk(const char *                   uplo,
     const char *                   trans,
     const dealii::types::blas_int *n,
     const dealii::types::blas_int *k,
     const float *                  alpha,
     const float *                  a,
     const dealii::types::blas_int *lda,
     const float *                  beta,
     float *                        c,
     const dealii::types::blas_int *ldc)
{
#ifdef DEAL_II_WITH_LAPACK
  ssyrk_(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
#else
  (void)uplo;
  (void)trans;
  (void)n;
  (void)k;
  (void)alpha;
  (void)a;
  (void)lda;
  (void)beta;
  (void)c;
  (void)ldc;
  Assert(false, LAPACKSupport::ExcMissing("ssyrk"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
syrk(const char *                   uplo,
     const char *                   trans,
     const dealii::types::blas_int *n,
     const dealii::types::blas_int *k,
     const double *                 alpha,
     const double *                 a,
     const dealii::types::blas_int *lda,
     const double *                 beta,
     double *                       c,
     const dealii::types::blas_int *ldc)
{
#ifdef DEAL_II_WITH_LAPACK
  dsyrk_(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
#else
  (void)uplo;
  (void)trans;
  (void)n;
  (void)k;
  (void)alpha;
  (void)a;
  (void)lda;
  (void)beta;
  (void)c;
  (void)ldc;
  Assert(false, LAPACKSupport::ExcMissing("dsyrk"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
syrk(const char *                   uplo,
     const char *                   trans,
     const dealii::types::blas_int *n,
     const dealii::types::blas_int *k,
     const std::complex<float> *    alpha,
     const std::complex<float> *    a,
     const dealii::types::blas_int *lda,
     const std::complex<float> *    beta,
     std::complex<float> *          c,
     const dealii::types::blas_int *ldc)
{
#ifdef DEAL_II_WITH_LAPACK
  csyrk_(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
#else
  (void)uplo;
  (void)trans;
  (void)n;
  (void)k;
  (void)alpha;
  (void)a;
  (void)lda;
  (void)beta;
  (void)c;
  (void)ldc;
  Assert(false, LAPACKSupport::ExcMissing("csyrk"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
syrk(const char *                   uplo,
     const char *                   trans,
     const dealii::types::blas_int *n,
     const dealii::types::blas_int *k,
     const std::complex<double> *   alpha,
     const std::complex<double> *   a,
     const dealii::types::blas_int *lda,
     const std::complex<double> *   beta,
     std::complex<double> *         c,
     const dealii::types::blas_int *ldc)
{
#ifdef DEAL_II_WITH_LAPACK
  zsyrk_(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
#else
  (void)uplo;
  (void)trans;
  (void)n;
  (void)k;
  (void)alpha;
  (void)a;
  (void)lda;
  (void)beta;
  (void)c;
  (void)ldc;
  Assert(false, LAPACKSupport::ExcMissing("zsyrk"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1, typename number2, typename number3>
inline void
trcon(const char *,
      const char *,
      const char *,
      const dealii::types::blas_int *,
      const number1 *,
      const dealii::types::blas_int *,
      number2 *,
      number3 *,
      dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
trcon(const char *                   norm,
      const char *                   uplo,
      const char *                   diag,
      const dealii::types::blas_int *n,
      const float *                  a,
      const dealii::types::blas_int *lda,
      float *                        rcond,
      float *                        work,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  strcon_(norm, uplo, diag, n, a, lda, rcond, work, iwork, info);
#else
  (void)norm;
  (void)uplo;
  (void)diag;
  (void)n;
  (void)a;
  (void)lda;
  (void)rcond;
  (void)work;
  (void)iwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("strcon"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
trcon(const char *                   norm,
      const char *                   uplo,
      const char *                   diag,
      const dealii::types::blas_int *n,
      const double *                 a,
      const dealii::types::blas_int *lda,
      double *                       rcond,
      double *                       work,
      dealii::types::blas_int *      iwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dtrcon_(norm, uplo, diag, n, a, lda, rcond, work, iwork, info);
#else
  (void)norm;
  (void)uplo;
  (void)diag;
  (void)n;
  (void)a;
  (void)lda;
  (void)rcond;
  (void)work;
  (void)iwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dtrcon"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
trcon(const char *                   norm,
      const char *                   uplo,
      const char *                   diag,
      const dealii::types::blas_int *n,
      const std::complex<float> *    a,
      const dealii::types::blas_int *lda,
      float *                        rcond,
      std::complex<float> *          work,
      float *                        rwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  ctrcon_(norm, uplo, diag, n, a, lda, rcond, work, rwork, info);
#else
  (void)norm;
  (void)uplo;
  (void)diag;
  (void)n;
  (void)a;
  (void)lda;
  (void)rcond;
  (void)work;
  (void)rwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("ctrcon"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
trcon(const char *                   norm,
      const char *                   uplo,
      const char *                   diag,
      const dealii::types::blas_int *n,
      const std::complex<double> *   a,
      const dealii::types::blas_int *lda,
      double *                       rcond,
      std::complex<double> *         work,
      double *                       rwork,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  ztrcon_(norm, uplo, diag, n, a, lda, rcond, work, rwork, info);
#else
  (void)norm;
  (void)uplo;
  (void)diag;
  (void)n;
  (void)a;
  (void)lda;
  (void)rcond;
  (void)work;
  (void)rwork;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("ztrcon"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1, typename number2>
inline void
trmv(const char *,
     const char *,
     const char *,
     const dealii::types::blas_int *,
     const number1 *,
     const dealii::types::blas_int *,
     number2 *,
     const dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
trmv(const char *                   uplo,
     const char *                   trans,
     const char *                   diag,
     const dealii::types::blas_int *n,
     const float *                  a,
     const dealii::types::blas_int *lda,
     float *                        x,
     const dealii::types::blas_int *incx)
{
#ifdef DEAL_II_WITH_LAPACK
  strmv_(uplo, trans, diag, n, a, lda, x, incx);
#else
  (void)uplo;
  (void)trans;
  (void)diag;
  (void)n;
  (void)a;
  (void)lda;
  (void)x;
  (void)incx;
  Assert(false, LAPACKSupport::ExcMissing("strmv"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
trmv(const char *                   uplo,
     const char *                   trans,
     const char *                   diag,
     const dealii::types::blas_int *n,
     const double *                 a,
     const dealii::types::blas_int *lda,
     double *                       x,
     const dealii::types::blas_int *incx)
{
#ifdef DEAL_II_WITH_LAPACK
  dtrmv_(uplo, trans, diag, n, a, lda, x, incx);
#else
  (void)uplo;
  (void)trans;
  (void)diag;
  (void)n;
  (void)a;
  (void)lda;
  (void)x;
  (void)incx;
  Assert(false, LAPACKSupport::ExcMissing("dtrmv"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
trmv(const char *                   uplo,
     const char *                   trans,
     const char *                   diag,
     const dealii::types::blas_int *n,
     const std::complex<float> *    a,
     const dealii::types::blas_int *lda,
     std::complex<float> *          x,
     const dealii::types::blas_int *incx)
{
#ifdef DEAL_II_WITH_LAPACK
  ctrmv_(uplo, trans, diag, n, a, lda, x, incx);
#else
  (void)uplo;
  (void)trans;
  (void)diag;
  (void)n;
  (void)a;
  (void)lda;
  (void)x;
  (void)incx;
  Assert(false, LAPACKSupport::ExcMissing("ctrmv"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
trmv(const char *                   uplo,
     const char *                   trans,
     const char *                   diag,
     const dealii::types::blas_int *n,
     const std::complex<double> *   a,
     const dealii::types::blas_int *lda,
     std::complex<double> *         x,
     const dealii::types::blas_int *incx)
{
#ifdef DEAL_II_WITH_LAPACK
  ztrmv_(uplo, trans, diag, n, a, lda, x, incx);
#else
  (void)uplo;
  (void)trans;
  (void)diag;
  (void)n;
  (void)a;
  (void)lda;
  (void)x;
  (void)incx;
  Assert(false, LAPACKSupport::ExcMissing("ztrmv"));
#endif // DEAL_II_WITH_LAPACK
}



template <typename number1, typename number2>
inline void
trtrs(const char *,
      const char *,
      const char *,
      const dealii::types::blas_int *,
      const dealii::types::blas_int *,
      const number1 *,
      const dealii::types::blas_int *,
      number2 *,
      const dealii::types::blas_int *,
      dealii::types::blas_int *)
{
  Assert(false, ExcNotImplemented());
}



inline void
trtrs(const char *                   uplo,
      const char *                   trans,
      const char *                   diag,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      const float *                  a,
      const dealii::types::blas_int *lda,
      float *                        b,
      const dealii::types::blas_int *ldb,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  strtrs_(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
#else
  (void)uplo;
  (void)trans;
  (void)diag;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("strtrs"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
trtrs(const char *                   uplo,
      const char *                   trans,
      const char *                   diag,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      const double *                 a,
      const dealii::types::blas_int *lda,
      double *                       b,
      const dealii::types::blas_int *ldb,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  dtrtrs_(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
#else
  (void)uplo;
  (void)trans;
  (void)diag;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("dtrtrs"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
trtrs(const char *                   uplo,
      const char *                   trans,
      const char *                   diag,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      const std::complex<float> *    a,
      const dealii::types::blas_int *lda,
      std::complex<float> *          b,
      const dealii::types::blas_int *ldb,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  ctrtrs_(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
#else
  (void)uplo;
  (void)trans;
  (void)diag;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("ctrtrs"));
#endif // DEAL_II_WITH_LAPACK
}



inline void
trtrs(const char *                   uplo,
      const char *                   trans,
      const char *                   diag,
      const dealii::types::blas_int *n,
      const dealii::types::blas_int *nrhs,
      const std::complex<double> *   a,
      const dealii::types::blas_int *lda,
      std::complex<double> *         b,
      const dealii::types::blas_int *ldb,
      dealii::types::blas_int *      info)
{
#ifdef DEAL_II_WITH_LAPACK
  ztrtrs_(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
#else
  (void)uplo;
  (void)trans;
  (void)diag;
  (void)n;
  (void)nrhs;
  (void)a;
  (void)lda;
  (void)b;
  (void)ldb;
  (void)info;
  Assert(false, LAPACKSupport::ExcMissing("ztrtrs"));
#endif // DEAL_II_WITH_LAPACK
}



// These overloads of lamch were added for 9.0 but are replaced by the
// automatically generated versions provided above. These are manual additions
// and should be kept even if this file is regenerated.
template <typename number>
DEAL_II_DEPRECATED inline void
lamch(const char *, number &)
{
  Assert(false, ExcNotImplemented());
}



DEAL_II_DEPRECATED
inline void
lamch(const char *cmach, float &precision)
{
#ifdef DEAL_II_WITH_LAPACK
  precision = slamch_(cmach);
#else
  (void)cmach;
  (void)precision;
  Assert(false, LAPACKSupport::ExcMissing("slamch"));
#endif // DEAL_II_WITH_LAPACK
}



DEAL_II_DEPRECATED
inline void
lamch(const char *cmach, double &precision)
{
#ifdef DEAL_II_WITH_LAPACK
  precision = dlamch_(cmach);
#else
  (void)cmach;
  (void)precision;
  Assert(false, LAPACKSupport::ExcMissing("dlamch"));
#endif // DEAL_II_WITH_LAPACK
}

DEAL_II_NAMESPACE_CLOSE

#endif
