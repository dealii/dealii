//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    This file was automatically generated from blas.h.in
//    See blastemplates in the deal.II contrib directory
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __LAPACK_TEMPLATES_H
#define __LAPACK_TEMPLATES_H

extern "C"
{
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
// Compute eigenvalues and vectors
void dgeev_ (char* jobvl, char* jobvr,
	     int* n, double* A, int* lda,
	     double* lambda_re, double* lambda_im,
	     double* vl, int* ldvl,
	     double* vr, int* ldva,
	     double* work, int* lwork,
	     int* info);
void sgeev_ (char* jobvl, char* jobvr,
	     int* n, float* A, int* lda,
	     float* lambda_re, float* lambda_im,
	     float* vl, int* ldvl,
	     float* vr, int* ldva,
	     float* work, int* lwork,
	     int* info);
// Compute eigenvalues and vectors (expert)
void dgeevx_ (char* balanc, char* jobvl, char* jobvr, char* sense,
	      int* n, double* A, int* lda,
	      double* lambda_re, double* lambda_im,
	      double* vl, int* ldvl,
	      double* vr, int* ldvr,
	      int* ilo, int* ihi,
	      double* scale, double* abnrm,
	      double* rconde, double* rcondv,
	      double* work, int* lwork,
	      int* iwork, int* info);
void sgeevx_ (char* balanc, char* jobvl, char* jobvr, char* sense,
	      int* n, float* A, int* lda,
	      float* lambda_re, float* lambda_im,
	      float* vl, int* ldvl,
	      float* vr, int* ldvr,
	      int* ilo, int* ihi,
	      float* scale, float* abnrm,
	      float* rconde, float* rcondv,
	      float* work, int* lwork,
	      int* iwork, int* info);
// Compute singular value decomposition
void dgesvd_ (int* jobu, int* jobvt,
	      int* n, int* m, double* A, int* lda,
	      double* s,
	      double* u, int* ldu,
	      double* vt, int* ldvt,
	      double* work, int* lwork,
	      int* info);
void sgesvd_ (int* jobu, int* jobvt,
	      int* n, int* m, float* A, int* lda,
	      float* s,
	      float* u, int* ldu,
	      float* vt, int* ldvt,
	      float* work, int* lwork,
	      int* info);

}



inline void
gemv (const char* trans, const int* m, const int* n, const double* alpha, const double* A, const int* lda, const double* x, const int* incx, const double* b, double* y, const int* incy)
{
  dgemv_ (trans,m,n,alpha,A,lda,x,incx,b,y,incy);
}


inline void
gemv (const char* trans, const int* m, const int* n, const float* alpha, const float* A, const int* lda, const float* x, const int* incx, const float* b, float* y, const int* incy)
{
  sgemv_ (trans,m,n,alpha,A,lda,x,incx,b,y,incy);
}


inline void
geev (char* jobvl, char* jobvr, int* n, double* A, int* lda, double* lambda_re, double* lambda_im, double* vl, int* ldvl, double* vr, int* ldva, double* work, int* lwork, int* info)
{
  dgeev_ (jobvl,jobvr,n,A,lda,lambda_re,lambda_im,vl,ldvl,vr,ldva,work,lwork,info);
}


inline void
geev (char* jobvl, char* jobvr, int* n, float* A, int* lda, float* lambda_re, float* lambda_im, float* vl, int* ldvl, float* vr, int* ldva, float* work, int* lwork, int* info)
{
  sgeev_ (jobvl,jobvr,n,A,lda,lambda_re,lambda_im,vl,ldvl,vr,ldva,work,lwork,info);
}


inline void
geevx (char* balanc, char* jobvl, char* jobvr, char* sense, int* n, double* A, int* lda, double* lambda_re, double* lambda_im, double* vl, int* ldvl, double* vr, int* ldvr, int* ilo, int* ihi, double* scale, double* abnrm, double* rconde, double* rcondv, double* work, int* lwork, int* iwork, int* info)
{
  dgeevx_ (balanc,jobvl,jobvr,sense,n,A,lda,lambda_re,lambda_im,vl,ldvl,vr,ldvr,ilo,ihi,scale,abnrm,rconde,rcondv,work,lwork,iwork,info);
}


inline void
geevx (char* balanc, char* jobvl, char* jobvr, char* sense, int* n, float* A, int* lda, float* lambda_re, float* lambda_im, float* vl, int* ldvl, float* vr, int* ldvr, int* ilo, int* ihi, float* scale, float* abnrm, float* rconde, float* rcondv, float* work, int* lwork, int* iwork, int* info)
{
  sgeevx_ (balanc,jobvl,jobvr,sense,n,A,lda,lambda_re,lambda_im,vl,ldvl,vr,ldvr,ilo,ihi,scale,abnrm,rconde,rcondv,work,lwork,iwork,info);
}


inline void
gesvd (int* jobu, int* jobvt, int* n, int* m, double* A, int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt, double* work, int* lwork, int* info)
{
  dgesvd_ (jobu,jobvt,n,m,A,lda,s,u,ldu,vt,ldvt,work,lwork,info);
}


inline void
gesvd (int* jobu, int* jobvt, int* n, int* m, float* A, int* lda, float* s, float* u, int* ldu, float* vt, int* ldvt, float* work, int* lwork, int* info)
{
  sgesvd_ (jobu,jobvt,n,m,A,lda,s,u,ldu,vt,ldvt,work,lwork,info);
}


#endif
