//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    This file was automatically generated from lapack_templates.h.in
//    See blastemplates in the deal.II contrib directory
//
//    Copyright (C) 2005, 2006, 2007, 2008, 2009 by the deal authors
//
//    This file is subject to QPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __LAPACK_TEMPLATES_H
#define __LAPACK_TEMPLATES_H

#include <base/config.h>
#include <lac/lapack_support.h>

extern "C"
{
// vector update of the form y += alpha*x with a scalar, x,y vectors
void daxpy_ (const int* n, const double* alpha, const double* x,
	     const int* incx, double* y, const int* incy);
void saxpy_ (const int* n, const float* alpha, const float* x,
	     const int* incx, float* y, const int* incy);
// General Matrix
// Matrix vector product
void dgemv_ (const char* trans, const int* m, const int* n,
	     const double* alpha, const double* A, const int* lda,
	     const double* x, const int* incx,
	     const double* b, double* y, const int* incy);
void sgemv_ (const char* trans, const int* m, const int* n,
	     const float* alpha, const float* A, const int* lda,
	     const float* x, const int* incx,
	     const float* b, float* y, const int* incy);
// Matrix matrix product
void dgemm_ (const char* transa, const char* transb,
	     const int* m, const int* n, const int* k,
	     const double* alpha, const double* A, const int* lda,
	     const double* B, const int* ldb,
	     const double* beta, double* C, const int* ldc);
void sgemm_ (const char* transa, const char* transb,
	     const int* m, const int* n, const int* k,
	     const float* alpha, const float* A, const int* lda,
	     const float* B, const int* ldb,
	     const float* beta, float* C, const int* ldc);
// Compute LU factorization
void dgetrf_ (const int* m, const int* n, double* A,
	      const int* lda, int* ipiv, int* info);
void sgetrf_ (const int* m, const int* n, float* A,
	      const int* lda, int* ipiv, int* info);
// Apply forward/backward substitution to LU factorization
void dgetrs_ (const char* trans, const int* n, const int* nrhs,
	      const double* A, const int* lda, const int* ipiv,
	      double* b, const int* ldb, int* info);
void sgetrs_ (const char* trans, const int* n, const int* nrhs,
	      const float* A, const int* lda, const int* ipiv,
	      float* b, const int* ldb, int* info);
// Invert matrix from LU factorization
void dgetri_ (const int* n, double* A, const int* lda, 
	      int* ipiv, double* inv_work, const int* lwork, int* info);
void sgetri_ (const int* n, float* A, const int* lda, 
	      int* ipiv, float* inv_work, const int* lwork, int* info);
// Compute QR factorization (Householder)
void dgeqrf_ (const int* m, const int* n, double* A,
	      const int* lda, double* tau, double* work,
	      const int* lwork, int* info);
void sgeqrf_ (const int* m, const int* n, float* A,
	      const int* lda, float* tau, float* work,
	      const int* lwork, int* info);
// Compute vector Q^T B, where Q is the result from dgeqrf_
void dormqr_ (const char* side, const char* trans, const int* m, 
	      const int* n, const int* k, const double* A, const int* lda, 
	      const double* tau, double* B, const int* ldb,
	      double* work, const int* lwork, int* info);
void sormqr_ (const char* side, const char* trans, const int* m, 
	      const int* n, const int* k, const float* A, const int* lda, 
	      const float* tau, float* B, const int* ldb,
	      float* work, const int* lwork, int* info);
// Compute matrix Q from the result of dgeqrf_
void dorgqr_ (const int* m, const int* n, const int* k, const double* A, 
	      const int* lda, const double* tau, double* work, const int* lwork, 
	      int* info);
void sorgqr_ (const int* m, const int* n, const int* k, const float* A, 
	      const int* lda, const float* tau, float* work, const int* lwork, 
	      int* info);
// Compute Rx = b
void dtrtrs_ (const char* uplo, const char* trans,
	      const char* diag, const int* n, const int* n_rhs,
	      const double* A, const int* lda, double* B, const int* ldb,
	      int* info);
void strtrs_ (const char* uplo, const char* trans,
	      const char* diag, const int* n, const int* n_rhs,
	      const float* A, const int* lda, float* B, const int* ldb,
	      int* info);
// Compute eigenvalues and vectors
void dgeev_ (const char* jobvl, const char* jobvr,
	     const int* n, double* A, const int* lda,
	     double* lambda_re, double* lambda_im,
	     double* vl, const int* ldvl,
	     double* vr, const int* ldva,
	     double* work, const int* lwork,
	     int* info);
void sgeev_ (const char* jobvl, const char* jobvr,
	     const int* n, float* A, const int* lda,
	     float* lambda_re, float* lambda_im,
	     float* vl, const int* ldvl,
	     float* vr, const int* ldva,
	     float* work, const int* lwork,
	     int* info);
// Compute eigenvalues and vectors (expert)
void dgeevx_ (const char* balanc, const char* jobvl, const char* jobvr,
	      const char* sense,
	      const int* n, double* A, const int* lda,
	      double* lambda_re, double* lambda_im,
	      double* vl, const int* ldvl,
	      double* vr, const int* ldvr,
	      int* ilo, int* ihi,
	      double* scale, double* abnrm,
	      double* rconde, double* rcondv,
	      double* work, const int* lwork,
	      int* iwork, int* info);
void sgeevx_ (const char* balanc, const char* jobvl, const char* jobvr,
	      const char* sense,
	      const int* n, float* A, const int* lda,
	      float* lambda_re, float* lambda_im,
	      float* vl, const int* ldvl,
	      float* vr, const int* ldvr,
	      int* ilo, int* ihi,
	      float* scale, float* abnrm,
	      float* rconde, float* rcondv,
	      float* work, const int* lwork,
	      int* iwork, int* info);
// Eigenvalues for a symmetric matrix
void dsyev_ (const char *jobz, const char *uplo, const int *n,
	     double *A, const int *lda, double *w,
	     double *work, const int *lwork, int *info);
void ssyev_ (const char *jobz, const char *uplo, const int *n,
	     float *A, const int *lda, float *w,
	     float *work, const int *lwork, int *info);
// Compute singular value decomposition
void dgesvd_ (int* jobu, int* jobvt,
	      const int* n, const int* m, double* A, const int* lda,
	      double* s,
	      double* u, const int* ldu,
	      double* vt, const int* ldvt,
	      double* work, const int* lwork,
	      int* info);
void sgesvd_ (int* jobu, int* jobvt,
	      const int* n, const int* m, float* A, const int* lda,
	      float* s,
	      float* u, const int* ldu,
	      float* vt, const int* ldvt,
	      float* work, const int* lwork,
	      int* info);
// Symmetric tridiagonal matrix
void dstev_ (const char* jobz, const int* n,
	     double* d, double* e, double* z,
	     const int* ldz, double* work,
	     int* info);
void sstev_ (const char* jobz, const int* n,
	     float* d, float* e, float* z,
	     const int* ldz, float* work,
	     int* info);

}


DEAL_II_NAMESPACE_OPEN

#ifdef HAVE_DAXPY_
inline void
axpy (const int* n, const double* alpha, const double* x, const int* incx, double* y, const int* incy)
{
  daxpy_ (n,alpha,x,incx,y,incy);
}
#else
inline void
axpy (const int*, const double*, const double*, const int*, double*, const int*)
{
  Assert (false, LAPACKSupport::ExcMissing("daxpy"));
}
#endif


#ifdef HAVE_SAXPY_
inline void
axpy (const int* n, const float* alpha, const float* x, const int* incx, float* y, const int* incy)
{
  saxpy_ (n,alpha,x,incx,y,incy);
}
#else
inline void
axpy (const int*, const float*, const float*, const int*, float*, const int*)
{
  Assert (false, LAPACKSupport::ExcMissing("saxpy"));
}
#endif


#ifdef HAVE_DGEMV_
inline void
gemv (const char* trans, const int* m, const int* n, const double* alpha, const double* A, const int* lda, const double* x, const int* incx, const double* b, double* y, const int* incy)
{
  dgemv_ (trans,m,n,alpha,A,lda,x,incx,b,y,incy);
}
#else
inline void
gemv (const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*)
{
  Assert (false, LAPACKSupport::ExcMissing("dgemv"));
}
#endif


#ifdef HAVE_SGEMV_
inline void
gemv (const char* trans, const int* m, const int* n, const float* alpha, const float* A, const int* lda, const float* x, const int* incx, const float* b, float* y, const int* incy)
{
  sgemv_ (trans,m,n,alpha,A,lda,x,incx,b,y,incy);
}
#else
inline void
gemv (const char*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*)
{
  Assert (false, LAPACKSupport::ExcMissing("sgemv"));
}
#endif


#ifdef HAVE_DGEMM_
inline void
gemm (const char* transa, const char* transb, const int* m, const int* n, const int* k, const double* alpha, const double* A, const int* lda, const double* B, const int* ldb, const double* beta, double* C, const int* ldc)
{
  dgemm_ (transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
}
#else
inline void
gemm (const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*)
{
  Assert (false, LAPACKSupport::ExcMissing("dgemm"));
}
#endif


#ifdef HAVE_SGEMM_
inline void
gemm (const char* transa, const char* transb, const int* m, const int* n, const int* k, const float* alpha, const float* A, const int* lda, const float* B, const int* ldb, const float* beta, float* C, const int* ldc)
{
  sgemm_ (transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
}
#else
inline void
gemm (const char*, const char*, const int*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*)
{
  Assert (false, LAPACKSupport::ExcMissing("sgemm"));
}
#endif


#ifdef HAVE_DGETRF_
inline void
getrf (const int* m, const int* n, double* A, const int* lda, int* ipiv, int* info)
{
  dgetrf_ (m,n,A,lda,ipiv,info);
}
#else
inline void
getrf (const int*, const int*, double*, const int*, int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("dgetrf"));
}
#endif


#ifdef HAVE_SGETRF_
inline void
getrf (const int* m, const int* n, float* A, const int* lda, int* ipiv, int* info)
{
  sgetrf_ (m,n,A,lda,ipiv,info);
}
#else
inline void
getrf (const int*, const int*, float*, const int*, int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("sgetrf"));
}
#endif


#ifdef HAVE_DGETRS_
inline void
getrs (const char* trans, const int* n, const int* nrhs, const double* A, const int* lda, const int* ipiv, double* b, const int* ldb, int* info)
{
  dgetrs_ (trans,n,nrhs,A,lda,ipiv,b,ldb,info);
}
#else
inline void
getrs (const char*, const int*, const int*, const double*, const int*, const int*, double*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("dgetrs"));
}
#endif


#ifdef HAVE_SGETRS_
inline void
getrs (const char* trans, const int* n, const int* nrhs, const float* A, const int* lda, const int* ipiv, float* b, const int* ldb, int* info)
{
  sgetrs_ (trans,n,nrhs,A,lda,ipiv,b,ldb,info);
}
#else
inline void
getrs (const char*, const int*, const int*, const float*, const int*, const int*, float*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("sgetrs"));
}
#endif


#ifdef HAVE_DGETRI_
inline void
getri (const int* n, double* A, const int* lda, int* ipiv, double* inv_work, const int* lwork, int* info)
{
  dgetri_ (n,A,lda,ipiv,inv_work,lwork,info);
}
#else
inline void
getri (const int*, double*, const int*, int*, double*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("dgetri"));
}
#endif


#ifdef HAVE_SGETRI_
inline void
getri (const int* n, float* A, const int* lda, int* ipiv, float* inv_work, const int* lwork, int* info)
{
  sgetri_ (n,A,lda,ipiv,inv_work,lwork,info);
}
#else
inline void
getri (const int*, float*, const int*, int*, float*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("sgetri"));
}
#endif


#ifdef HAVE_DGEQRF_
inline void
geqrf (const int* m, const int* n, double* A, const int* lda, double* tau, double* work, const int* lwork, int* info)
{
  dgeqrf_ (m,n,A,lda,tau,work,lwork,info);
}
#else
inline void
geqrf (const int*, const int*, double*, const int*, double*, double*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("dgeqrf"));
}
#endif


#ifdef HAVE_SGEQRF_
inline void
geqrf (const int* m, const int* n, float* A, const int* lda, float* tau, float* work, const int* lwork, int* info)
{
  sgeqrf_ (m,n,A,lda,tau,work,lwork,info);
}
#else
inline void
geqrf (const int*, const int*, float*, const int*, float*, float*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("sgeqrf"));
}
#endif


#ifdef HAVE_DORMQR_
inline void
ormqr (const char* side, const char* trans, const int* m, const int* n, const int* k, const double* A, const int* lda, const double* tau, double* B, const int* ldb, double* work, const int* lwork, int* info)
{
  dormqr_ (side,trans,m,n,k,A,lda,tau,B,ldb,work,lwork,info);
}
#else
inline void
ormqr (const char*, const char*, const int*, const int*, const int*, const double*, const int*, const double*, double*, const int*, double*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("dormqr"));
}
#endif


#ifdef HAVE_SORMQR_
inline void
ormqr (const char* side, const char* trans, const int* m, const int* n, const int* k, const float* A, const int* lda, const float* tau, float* B, const int* ldb, float* work, const int* lwork, int* info)
{
  sormqr_ (side,trans,m,n,k,A,lda,tau,B,ldb,work,lwork,info);
}
#else
inline void
ormqr (const char*, const char*, const int*, const int*, const int*, const float*, const int*, const float*, float*, const int*, float*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("sormqr"));
}
#endif


#ifdef HAVE_DORGQR_
inline void
orgqr (const int* m, const int* n, const int* k, const double* A, const int* lda, const double* tau, double* work, const int* lwork, int* info)
{
  dorgqr_ (m,n,k,A,lda,tau,work,lwork,info);
}
#else
inline void
orgqr (const int*, const int*, const int*, const double*, const int*, const double*, double*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("dorgqr"));
}
#endif


#ifdef HAVE_SORGQR_
inline void
orgqr (const int* m, const int* n, const int* k, const float* A, const int* lda, const float* tau, float* work, const int* lwork, int* info)
{
  sorgqr_ (m,n,k,A,lda,tau,work,lwork,info);
}
#else
inline void
orgqr (const int*, const int*, const int*, const float*, const int*, const float*, float*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("sorgqr"));
}
#endif


#ifdef HAVE_DTRTRS_
inline void
trtrs (const char* uplo, const char* trans, const char* diag, const int* n, const int* n_rhs, const double* A, const int* lda, double* B, const int* ldb, int* info)
{
  dtrtrs_ (uplo,trans,diag,n,n_rhs,A,lda,B,ldb,info);
}
#else
inline void
trtrs (const char*, const char*, const char*, const int*, const int*, const double*, const int*, double*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("dtrtrs"));
}
#endif


#ifdef HAVE_STRTRS_
inline void
trtrs (const char* uplo, const char* trans, const char* diag, const int* n, const int* n_rhs, const float* A, const int* lda, float* B, const int* ldb, int* info)
{
  strtrs_ (uplo,trans,diag,n,n_rhs,A,lda,B,ldb,info);
}
#else
inline void
trtrs (const char*, const char*, const char*, const int*, const int*, const float*, const int*, float*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("strtrs"));
}
#endif


#ifdef HAVE_DGEEV_
inline void
geev (const char* jobvl, const char* jobvr, const int* n, double* A, const int* lda, double* lambda_re, double* lambda_im, double* vl, const int* ldvl, double* vr, const int* ldva, double* work, const int* lwork, int* info)
{
  dgeev_ (jobvl,jobvr,n,A,lda,lambda_re,lambda_im,vl,ldvl,vr,ldva,work,lwork,info);
}
#else
inline void
geev (const char*, const char*, const int*, double*, const int*, double*, double*, double*, const int*, double*, const int*, double*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("dgeev"));
}
#endif


#ifdef HAVE_SGEEV_
inline void
geev (const char* jobvl, const char* jobvr, const int* n, float* A, const int* lda, float* lambda_re, float* lambda_im, float* vl, const int* ldvl, float* vr, const int* ldva, float* work, const int* lwork, int* info)
{
  sgeev_ (jobvl,jobvr,n,A,lda,lambda_re,lambda_im,vl,ldvl,vr,ldva,work,lwork,info);
}
#else
inline void
geev (const char*, const char*, const int*, float*, const int*, float*, float*, float*, const int*, float*, const int*, float*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("sgeev"));
}
#endif


#ifdef HAVE_DGEEVX_
inline void
geevx (const char* balanc, const char* jobvl, const char* jobvr, const char* sense, const int* n, double* A, const int* lda, double* lambda_re, double* lambda_im, double* vl, const int* ldvl, double* vr, const int* ldvr, int* ilo, int* ihi, double* scale, double* abnrm, double* rconde, double* rcondv, double* work, const int* lwork, int* iwork, int* info)
{
  dgeevx_ (balanc,jobvl,jobvr,sense,n,A,lda,lambda_re,lambda_im,vl,ldvl,vr,ldvr,ilo,ihi,scale,abnrm,rconde,rcondv,work,lwork,iwork,info);
}
#else
inline void
geevx (const char*, const char*, const char*, const char*, const int*, double*, const int*, double*, double*, double*, const int*, double*, const int*, int*, int*, double*, double*, double*, double*, double*, const int*, int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("dgeevx"));
}
#endif


#ifdef HAVE_SGEEVX_
inline void
geevx (const char* balanc, const char* jobvl, const char* jobvr, const char* sense, const int* n, float* A, const int* lda, float* lambda_re, float* lambda_im, float* vl, const int* ldvl, float* vr, const int* ldvr, int* ilo, int* ihi, float* scale, float* abnrm, float* rconde, float* rcondv, float* work, const int* lwork, int* iwork, int* info)
{
  sgeevx_ (balanc,jobvl,jobvr,sense,n,A,lda,lambda_re,lambda_im,vl,ldvl,vr,ldvr,ilo,ihi,scale,abnrm,rconde,rcondv,work,lwork,iwork,info);
}
#else
inline void
geevx (const char*, const char*, const char*, const char*, const int*, float*, const int*, float*, float*, float*, const int*, float*, const int*, int*, int*, float*, float*, float*, float*, float*, const int*, int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("sgeevx"));
}
#endif


#ifdef HAVE_DSYEV_
inline void
syev (const char *jobz, const char *uplo, const int *n, double *A, const int *lda, double *w, double *work, const int *lwork, int *info)
{
  dsyev_ (char*jobz,char*uplo,int*n,double*A,int*lda,double*w,double*work,int*lwork,int*info);
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
  ssyev_ (char*jobz,char*uplo,int*n,double*A,int*lda,double*w,double*work,int*lwork,int*info);
}
#else
inline void
syev (const char *, const char *, const int *, float *, const int *, float *, float *, const int *, int *)
{
  Assert (false, LAPACKSupport::ExcMissing("ssyev"));
}
#endif


#ifdef HAVE_DGESVD_
inline void
gesvd (int* jobu, int* jobvt, const int* n, const int* m, double* A, const int* lda, double* s, double* u, const int* ldu, double* vt, const int* ldvt, double* work, const int* lwork, int* info)
{
  dgesvd_ (jobu,jobvt,n,m,A,lda,s,u,ldu,vt,ldvt,work,lwork,info);
}
#else
inline void
gesvd (int*, int*, const int*, const int*, double*, const int*, double*, double*, const int*, double*, const int*, double*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("dgesvd"));
}
#endif


#ifdef HAVE_SGESVD_
inline void
gesvd (int* jobu, int* jobvt, const int* n, const int* m, float* A, const int* lda, float* s, float* u, const int* ldu, float* vt, const int* ldvt, float* work, const int* lwork, int* info)
{
  sgesvd_ (jobu,jobvt,n,m,A,lda,s,u,ldu,vt,ldvt,work,lwork,info);
}
#else
inline void
gesvd (int*, int*, const int*, const int*, float*, const int*, float*, float*, const int*, float*, const int*, float*, const int*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("sgesvd"));
}
#endif


#ifdef HAVE_DSTEV_
inline void
stev (const char* jobz, const int* n, double* d, double* e, double* z, const int* ldz, double* work, int* info)
{
  dstev_ (jobz,n,d,e,z,ldz,work,info);
}
#else
inline void
stev (const char*, const int*, double*, double*, double*, const int*, double*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("dstev"));
}
#endif


#ifdef HAVE_SSTEV_
inline void
stev (const char* jobz, const int* n, float* d, float* e, float* z, const int* ldz, float* work, int* info)
{
  sstev_ (jobz,n,d,e,z,ldz,work,info);
}
#else
inline void
stev (const char*, const int*, float*, float*, float*, const int*, float*, int*)
{
  Assert (false, LAPACKSupport::ExcMissing("sstev"));
}
#endif

DEAL_II_NAMESPACE_CLOSE

#endif
