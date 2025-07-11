// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_scalapack_templates_h
#define dealii_scalapack_templates_h


#include <deal.II/base/config.h>

// This file does not actually import anything into namespace dealii,
// but to avoid it being completely empty to some of our scripts, we
// need to make sure it opens and closes the namespace at least once.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE // Do not convert for module purposes


#ifdef DEAL_II_WITH_SCALAPACK

#  include <deal.II/base/mpi.h>
#  include <deal.II/base/mpi.templates.h>

#  ifdef DEAL_II_HAVE_FP_EXCEPTIONS
#    include <cfenv>
#  endif


  // useful examples:
  // https://stackoverflow.com/questions/14147705/cholesky-decomposition-scalapack-error/14203864
  // http://icl.cs.utk.edu/lapack-forum/viewtopic.php?t=139   // second post by
  // Julien Langou
  // https://andyspiros.wordpress.com/2011/07/08/an-example-of-blacs-with-c/
  // http://qboxcode.org/trac/browser/qb/tags/rel1_63_4/src/Matrix.C
  // https://gitlab.phys.ethz.ch/lwossnig/lecture/blob/a534f562dfb2ad5c564abe5c2356d5d956fb7218/examples/mpi/scalapack.cpp
  // https://github.com/elemental/Elemental/blob/master/src/core/imports/scalapack.cpp
  // https://scicomp.stackexchange.com/questions/7766/performance-optimization-or-tuning-possible-for-scalapack-gemm
  //
  // info:
  // http://www.netlib.org/scalapack/slug/index.html       // User guide
  // http://www.netlib.org/scalapack/slug/node135.html // How to Measure Errors

  extern "C"
{
  /* Basic Linear Algebra Communication Subprograms (BLACS) declarations */
  // https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dinitb.htm#dinitb

  /**
   * Determine how many processes are available and the current process rank.
   *
   * https://www.ibm.com/support/knowledgecenter/en/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dbpnf.htm
   */
  void Cblacs_pinfo(int *rank, int *nprocs);

  /**
   * Return internal BLACS value in @p val based on the input @p what and @p icontxt.
   * The most common use is in retrieving a default system context (@p what = 0, @p icontxt is ignored)
   * to be used in BLACS_GRIDINIT or BLACS_GRIDMAP.
   *
   * https://www.ibm.com/support/knowledgecenter/en/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dbget.htm
   */
  void Cblacs_get(int icontxt, int what, int *val);

  /**
   * Map the processes sequentially in row-major or column-major order
   * into the process grid. Input arguments must be the same on every process.
   *
   * On return, @p context is the integer handle to the BLACS context,
   * whereas on entry it is a system context to be used in creating the
   * BLACS context.
   *
   * https://www.ibm.com/support/knowledgecenter/en/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dbint.htm
   */
  void Cblacs_gridinit(int        *context,
                       const char *order,
                       int         grid_height,
                       int         grid_width);

  /**
   * Return the process row and column index.
   *
   * https://www.ibm.com/support/knowledgecenter/en/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dbinfo.htm
   */
  void Cblacs_gridinfo(int  context,
                       int *grid_height,
                       int *grid_width,
                       int *grid_row,
                       int *grid_col);

  /**
   * Given the system process number, return the row and column coordinates in
   * the BLACS' process grid.
   */
  void Cblacs_pcoord(int ictxt, int pnum, int *prow, int *pcol);

  /**
   * Release a BLACS context.
   */
  void Cblacs_gridexit(int context);

  /**
   * This routines holds up execution of all processes within the indicated
   * scope until they have all called the routine.
   */
  void Cblacs_barrier(int, const char *);

  /**
   * Free all BLACS contexts and releases all allocated memory.
   */
  void Cblacs_exit(int error_code);

  /**
   * Receives a message from a process @prsrc, @p csrc into a general rectangular matrix.
   *
   * https://software.intel.com/en-us/mkl-developer-reference-c-gerv2d
   */
  void Cdgerv2d(
    int context, int M, int N, double *A, int lda, int rsrc, int csrc);
  void Csgerv2d(
    int context, int M, int N, float *A, int lda, int rsrc, int csrc);

  /**
   * Sends the general rectangular matrix A to the destination
   * process @p rdest @p cdest in the process grid.
   *
   * https://software.intel.com/en-us/mkl-developer-reference-c-2018-beta-gesd2d
   */
  void Cdgesd2d(
    int context, int M, int N, double *A, int lda, int rdest, int cdest);
  void Csgesd2d(
    int context, int M, int N, float *A, int lda, int rdest, int cdest);

  /**
   * Get BLACS context from MPI @p comm.
   */
  int Csys2blacs_handle(MPI_Comm comm);

  /**
   * Compute how many rows and columns each process owns (NUMber of Rows Or
   * Columns).
   *
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dnumy.htm
   */
  int numroc_(const int *n,
              const int *nb,
              const int *iproc,
              const int *isproc,
              const int *nprocs);

  /**
   * Compute the Cholesky factorization of an N-by-N real
   * symmetric positive definite distributed matrix sub( A ) denoting
   * A(IA:IA+N-1, JA:JA+N-1).
   *
   * http://www.netlib.org/scalapack/explore-html/d5/d9e/pdpotrf_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lpotrf.htm
   */
  void pdpotrf_(const char *UPLO,
                const int  *N,
                double     *A,
                const int  *IA,
                const int  *JA,
                const int  *DESCA,
                int        *INFO);
  void pspotrf_(const char *UPLO,
                const int  *N,
                float      *A,
                const int  *IA,
                const int  *JA,
                const int  *DESCA,
                int        *INFO);

  /**
   * Computes an LU factorization of a general distributed matrix sub( A )
   * using partial pivoting with row interchanges.
   *
   * http://www.netlib.org/scalapack/explore-html/df/dfe/pdgetrf_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/en/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lgetrf.htm
   */
  void pdgetrf_(const int *m,
                const int *n,
                double    *A,
                const int *IA,
                const int *JA,
                const int *DESCA,
                int       *ipiv,
                int       *INFO);
  void psgetrf_(const int *m,
                const int *n,
                float     *A,
                const int *IA,
                const int *JA,
                const int *DESCA,
                int       *ipiv,
                int       *INFO);

  /**
   * Compute the inverse of a real symmetric positive definite
   * distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1) using the
   * Cholesky factorization sub( A ) = U**T*U or L*L**T computed by
   * PDPOTRF.
   *
   * http://www.netlib.org/scalapack/explore-html/d2/d44/pdpotri_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lpotri.htm
   * https://software.intel.com/en-us/mkl-developer-reference-c-p-potri
   */
  void pdpotri_(const char *UPLO,
                const int  *N,
                double     *A,
                const int  *IA,
                const int  *JA,
                const int  *DESCA,
                int        *INFO);
  void pspotri_(const char *UPLO,
                const int  *N,
                float      *A,
                const int  *IA,
                const int  *JA,
                const int  *DESCA,
                int        *INFO);

  /**
   * PDGETRI computes the inverse of a distributed matrix using the LU
   * factorization computed by PDGETRF. This method inverts U and then
   * computes the inverse of sub( A ) = A(IA:IA+N-1,JA:JA+N-1) denoted
   * InvA by solving the system InvA*L = inv(U) for InvA.
   *
   * http://www.netlib.org/scalapack/explore-html/d3/df3/pdgetri_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lgetri.htm
   */
  void pdgetri_(const int *N,
                double    *A,
                const int *IA,
                const int *JA,
                const int *DESCA,
                const int *ipiv,
                double    *work,
                int       *lwork,
                int       *iwork,
                int       *liwork,
                int       *info);
  void psgetri_(const int *N,
                float     *A,
                const int *IA,
                const int *JA,
                const int *DESCA,
                const int *ipiv,
                float     *work,
                int       *lwork,
                int       *iwork,
                int       *liwork,
                int       *info);


  /**
   * PDTRTRI computes the inverse of a upper or lower triangular
   * distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1).
   *
   * http://www.netlib.org/scalapack/explore-html/d9/dc0/pdtrtri_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lpdtri.htm
   * https://software.intel.com/en-us/mkl-developer-reference-c-p-trtri
   */
  void pdtrtri_(const char *UPLO,
                const char *DIAG,
                const int  *N,
                double     *A,
                const int  *IA,
                const int  *JA,
                const int  *DESCA,
                int        *INFO);
  void pstrtri_(const char *UPLO,
                const char *DIAG,
                const int  *N,
                float      *A,
                const int  *IA,
                const int  *JA,
                const int  *DESCA,
                int        *INFO);

  /**
   * Estimate the reciprocal of the condition number (in the
   * l1-norm) of a real symmetric positive definite distributed matrix
   * using the Cholesky factorization.
   *
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lpocon.htm#lpocon
   * http://www.netlib.org/scalapack/explore-html/d4/df7/pdpocon_8f.html
   * https://software.intel.com/en-us/mkl-developer-reference-fortran-pocon
   */
  void pdpocon_(const char   *uplo,
                const int    *N,
                const double *A,
                const int    *IA,
                const int    *JA,
                const int    *DESCA,
                const double *ANORM,
                double       *RCOND,
                double       *WORK,
                const int    *LWORK,
                int          *IWORK,
                const int    *LIWORK,
                int          *INFO);
  void pspocon_(const char  *uplo,
                const int   *N,
                const float *A,
                const int   *IA,
                const int   *JA,
                const int   *DESCA,
                const float *ANORM,
                float       *RCOND,
                float       *WORK,
                const int   *LWORK,
                int         *IWORK,
                const int   *LIWORK,
                int         *INFO);

  /**
   * Norm of a real symmetric matrix
   *
   * http://www.netlib.org/scalapack/explore-html/dd/d12/pdlansy_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_pdlansy.htm#pdlansy
   */
  double pdlansy_(const char   *norm,
                  const char   *uplo,
                  const int    *N,
                  const double *A,
                  const int    *IA,
                  const int    *JA,
                  const int    *DESCA,
                  double       *work);
  float  pslansy_(const char  *norm,
                  const char  *uplo,
                  const int   *N,
                  const float *A,
                  const int   *IA,
                  const int   *JA,
                  const int   *DESCA,
                  float       *work);

  /**
   * Compute the Least Common Multiple (LCM) of two positive integers @p M and @p N.
   * In fact the routine Compute the greatest common divisor (GCD) and
   * use the fact that M*N = GCD*LCM.
   *
   * http://www.netlib.org/scalapack/explore-html/d0/d9b/ilcm_8f_source.html
   */
  int ilcm_(const int *M, const int *N);

  /**
   * Return the ceiling of the division of two integers.
   *
   * http://www.netlib.org/scalapack/explore-html/df/d07/iceil_8f_source.html
   */
  int iceil_(const int *i1, const int *i2);

  /**
   * Initialize the descriptor vector with the 8 input arguments
   */
  void descinit_(int       *desc,
                 const int *m,
                 const int *n,
                 const int *mb,
                 const int *nb,
                 const int *irsrc,
                 const int *icsrc,
                 const int *ictxt,
                 const int *lld,
                 int       *info);

  /**
   * Compute the global index of a distributed matrix entry
   * pointed to by the local index @p indxloc of the process indicated by
   * @p iproc.
   *
   * @param indxloc The local index of the distributed matrix entry.
   * @param nb Block size, size of the blocks the distributed matrix is split
   * into.
   * @param iproc The coordinate of the process whose local array row or column
   * is to be determined
   * @param isrcproc  The coordinate of the process that possesses the first
   * row/column of the distributed matrix
   * @param nprocs The total number processes over which the distributed matrix
   * is distributed
   */
  int indxl2g_(const int *indxloc,
               const int *nb,
               const int *iproc,
               const int *isrcproc,
               const int *nprocs);

  /**
   * Compute the solution to a real system of linear equations
   */
  void pdgesv_(const int *n,
               const int *nrhs,
               double    *A,
               const int *ia,
               const int *ja,
               const int *desca,
               int       *ipiv,
               double    *B,
               const int *ib,
               const int *jb,
               const int *descb,
               int       *info);
  void psgesv_(const int *n,
               const int *nrhs,
               float     *A,
               const int *ia,
               const int *ja,
               const int *desca,
               int       *ipiv,
               float     *B,
               const int *ib,
               const int *jb,
               const int *descb,
               int       *info);

  /**
   * Perform one of the matrix-matrix operations:
   * @f{align*}{
   * \mathrm{sub}(C) &\dealcoloneq \alpha op(\mathrm{sub}(A))op(\mathrm{sub}(B))
   *                            + \beta \mathrm{sub}(C), \\
   * \mathrm{sub}(C) &\dealcoloneq \alpha op(\mathrm{sub}(A))op(\mathrm{sub}(B))
   *                            + beta sub(C),
   * @f}
   * where
   * $\mathrm{sub}(C)$ denotes C(IC:IC+M-1,JC:JC+N-1),  and, $op(X)$ is one of
   * $op(X) = X$ or $op(X) = X^T$.
   */
  void pdgemm_(const char   *transa,
               const char   *transb,
               const int    *m,
               const int    *n,
               const int    *k,
               const double *alpha,
               const double *A,
               const int    *IA,
               const int    *JA,
               const int    *DESCA,
               const double *B,
               const int    *IB,
               const int    *JB,
               const int    *DESCB,
               const double *beta,
               double       *C,
               const int    *IC,
               const int    *JC,
               const int    *DESCC);
  void psgemm_(const char  *transa,
               const char  *transb,
               const int   *m,
               const int   *n,
               const int   *k,
               const float *alpha,
               const float *A,
               const int   *IA,
               const int   *JA,
               const int   *DESCA,
               const float *B,
               const int   *IB,
               const int   *JB,
               const int   *DESCB,
               const float *beta,
               float       *C,
               const int   *IC,
               const int   *JC,
               const int   *DESCC);

  /**
   * Return the value of the one norm, or the Frobenius norm, or the infinity
   * norm, or the element of largest absolute value of a distributed matrix
   */
  double pdlange_(const char   *norm,
                  const int    *m,
                  const int    *n,
                  const double *A,
                  const int    *ia,
                  const int    *ja,
                  const int    *desca,
                  double       *work);
  float  pslange_(const char  *norm,
                  const int   *m,
                  const int   *n,
                  const float *A,
                  const int   *ia,
                  const int   *ja,
                  const int   *desca,
                  float       *work);

  /**
   * Compute the process coordinate which possesses the entry of a
   * distributed matrix specified by a global index
   */
  int indxg2p_(const int *glob,
               const int *nb,
               const int *iproc,
               const int *isproc,
               const int *nprocs);

  /**
   * Compute all eigenvalues and, optionally, eigenvectors of a real symmetric
   * matrix A by calling the recommended sequence of ScaLAPACK routines. In its
   * present form, the routine assumes a homogeneous system and makes no checks
   * for consistency of the eigenvalues or eigenvectors across the different
   * processes. Because of this, it is possible that a heterogeneous system may
   * return incorrect results without any error messages.
   *
   * http://www.netlib.org/scalapack/explore-html/d0/d1a/pdsyev_8f.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lsyev.htm#lsyev
   */
  void pdsyev_(const char *jobz,
               const char *uplo,
               const int  *m,
               double     *A,
               const int  *ia,
               const int  *ja,
               int        *desca,
               double     *w,
               double     *z,
               const int  *iz,
               const int  *jz,
               int        *descz,
               double     *work,
               const int  *lwork,
               int        *info);
  void pssyev_(const char *jobz,
               const char *uplo,
               const int  *m,
               float      *A,
               const int  *ia,
               const int  *ja,
               int        *desca,
               float      *w,
               float      *z,
               const int  *iz,
               const int  *jz,
               int        *descz,
               float      *work,
               const int  *lwork,
               int        *info);

  /**
   * Copy all or a part of a distributed matrix A to another distributed matrix
   * B. No communication is performed, pdlacpy performs a local copy
   * $\mathrm{sub}(A) \dealcoloneq \mathrm{sub}(B)$, where $\mathrm{sub}(A)$
   * denotes $A(ia:ia+m-1, ja:ja+n-1)$ and $\mathrm{sub}(B)$ denotes
   * $B(ib:ib+m-1, jb:jb+n-1)$.
   */
  void pdlacpy_(const char   *uplo,
                const int    *m,
                const int    *n,
                const double *A,
                const int    *ia,
                const int    *ja,
                const int    *desca,
                double       *B,
                const int    *ib,
                const int    *jb,
                const int    *descb);
  void pslacpy_(const char  *uplo,
                const int   *m,
                const int   *n,
                const float *A,
                const int   *ia,
                const int   *ja,
                const int   *desca,
                float       *B,
                const int   *ib,
                const int   *jb,
                const int   *descb);

  /**
   * Copies the content of a general rectangular distributed matrix @p A to another distributed matrix @p B
   * It is not required that the matrices A and B have the same process grid or
   * block size, e.g. copying a matrix from a one-dimensional to a
   * two-dimensional process grid
   * @p ictxt is a context which is at least a union of all processes in context
   * A and B
   */
  void pdgemr2d_(const int    *m,
                 const int    *n,
                 const double *A,
                 const int    *ia,
                 const int    *ja,
                 const int    *desca,
                 double       *B,
                 const int    *ib,
                 const int    *jb,
                 const int    *descb,
                 const int    *ictxt);
  void psgemr2d_(const int   *m,
                 const int   *n,
                 const float *A,
                 const int   *ia,
                 const int   *ja,
                 const int   *desca,
                 float       *B,
                 const int   *ib,
                 const int   *jb,
                 const int   *descb,
                 const int   *ictxt);

  /**
   * helper routines determining machine precision
   */
  double pdlamch_(const int *ictxt, const char *cmach);
  float  pslamch_(const int *ictxt, const char *cmach);


  /**
   * psyevx computes selected eigenvalues and, optionally, eigenvectors
   * of a real symmetric matrix A. Eigenvalues/vectors can be selected by
   * specifying a range of values or a range of indices for the desired
   * eigenvalues.
   */
  void pdsyevx_(const char   *jobz,
                const char   *range,
                const char   *uplo,
                const int    *n,
                double       *A,
                const int    *ia,
                const int    *ja,
                const int    *desca,
                const double *VL,
                const double *VU,
                const int    *il,
                const int    *iu,
                const double *abstol,
                const int    *m,
                const int    *nz,
                double       *w,
                double       *orfac,
                double       *Z,
                const int    *iz,
                const int    *jz,
                const int    *descz,
                double       *work,
                int          *lwork,
                int          *iwork,
                int          *liwork,
                int          *ifail,
                int          *iclustr,
                double       *gap,
                int          *info);
  void pssyevx_(const char  *jobz,
                const char  *range,
                const char  *uplo,
                const int   *n,
                float       *A,
                const int   *ia,
                const int   *ja,
                const int   *desca,
                const float *VL,
                const float *VU,
                const int   *il,
                const int   *iu,
                const float *abstol,
                const int   *m,
                const int   *nz,
                float       *w,
                float       *orfac,
                float       *Z,
                const int   *iz,
                const int   *jz,
                const int   *descz,
                float       *work,
                int         *lwork,
                int         *iwork,
                int         *liwork,
                int         *ifail,
                int         *iclustr,
                float       *gap,
                int         *info);

  /*
   * PDGESVD computes the singular value decomposition (SVD) of an
   * M-by-N matrix A, optionally computing the left and/or right
   * singular vectors
   */
  void pdgesvd_(const char *jobu,
                const char *jobvt,
                const int  *m,
                const int  *n,
                double     *A,
                const int  *ia,
                const int  *ja,
                const int  *desca,
                double     *S,
                double     *U,
                const int  *iu,
                const int  *ju,
                const int  *descu,
                double     *VT,
                const int  *ivt,
                const int  *jvt,
                const int  *descvt,
                double     *work,
                int        *lwork,
                int        *info);
  void psgesvd_(const char *jobu,
                const char *jobvt,
                const int  *m,
                const int  *n,
                float      *A,
                const int  *ia,
                const int  *ja,
                const int  *desca,
                float      *S,
                float      *U,
                const int  *iu,
                const int  *ju,
                const int  *descu,
                float      *VT,
                const int  *ivt,
                const int  *jvt,
                const int  *descvt,
                float      *work,
                int        *lwork,
                int        *info);

  /*
   * P_GELS solves overdetermined or underdetermined real linear
   * systems involving an M-by-N matrix A, or its transpose,
   * using a QR or LQ factorization of A.  It is assumed that A has full rank.
   */
  void pdgels_(const char *trans,
               const int  *m,
               const int  *n,
               const int  *nrhs,
               double     *A,
               const int  *ia,
               const int  *ja,
               const int  *desca,
               double     *B,
               const int  *ib,
               const int  *jb,
               const int  *descb,
               double     *work,
               int        *lwork,
               int        *info);
  void psgels_(const char *trans,
               const int  *m,
               const int  *n,
               const int  *nrhs,
               float      *A,
               const int  *ia,
               const int  *ja,
               const int  *desca,
               float      *B,
               const int  *ib,
               const int  *jb,
               const int  *descb,
               float      *work,
               int        *lwork,
               int        *info);

  /*
   * Perform matrix sum:
   * @f{equation*}{
   * C \dealcoloneq \beta C + \alpha op(A),
   * @f
   * where $op(A)$ denotes either $op(A) = A$ or $op(A)=A^T$.
   */
  void pdgeadd_(const char   *transa,
                const int    *m,
                const int    *n,
                const double *alpha,
                const double *A,
                const int    *IA,
                const int    *JA,
                const int    *DESCA,
                const double *beta,
                double       *C,
                const int    *IC,
                const int    *JC,
                const int    *DESCC);
  void psgeadd_(const char  *transa,
                const int   *m,
                const int   *n,
                const float *alpha,
                const float *A,
                const int   *IA,
                const int   *JA,
                const int   *DESCA,
                const float *beta,
                float       *C,
                const int   *IC,
                const int   *JC,
                const int   *DESCC);

  /**
   * Routine to transpose a matrix:
   * C = beta C + alpha A^T
   */
  void pdtran_(const int    *m,
               const int    *n,
               const double *alpha,
               const double *A,
               const int    *IA,
               const int    *JA,
               const int    *DESCA,
               const double *beta,
               double       *C,
               const int    *IC,
               const int    *JC,
               const int    *DESCC);
  void pstran_(const int   *m,
               const int   *n,
               const float *alpha,
               const float *A,
               const int   *IA,
               const int   *JA,
               const int   *DESCA,
               const float *beta,
               float       *C,
               const int   *IC,
               const int   *JC,
               const int   *DESCC);

  /**
   * psyevr computes selected eigenvalues and, optionally, eigenvectors
   * of a real symmetric matrix A using a parallel implementation of the MRR
   * algorithm. Eigenvalues/vectors can be selected by specifying a range of
   * values or a range of indices for the desired eigenvalues.
   */
  void pdsyevr_(const char   *jobz,
                const char   *range,
                const char   *uplo,
                const int    *n,
                double       *A,
                const int    *IA,
                const int    *JA,
                const int    *DESCA,
                const double *VL,
                const double *VU,
                const int    *IL,
                const int    *IU,
                int          *m,
                int          *nz,
                double       *w,
                double       *Z,
                const int    *IZ,
                const int    *JZ,
                const int    *DESCZ,
                double       *work,
                int          *lwork,
                int          *iwork,
                int          *liwork,
                int          *info);
  void pssyevr_(const char  *jobz,
                const char  *range,
                const char  *uplo,
                const int   *n,
                float       *A,
                const int   *IA,
                const int   *JA,
                const int   *DESCA,
                const float *VL,
                const float *VU,
                const int   *IL,
                const int   *IU,
                int         *m,
                int         *nz,
                float       *w,
                float       *Z,
                const int   *IZ,
                const int   *JZ,
                const int   *DESCZ,
                float       *work,
                int         *lwork,
                int         *iwork,
                int         *liwork,
                int         *info);
}



/*
 * In the following we have template wrappers for the ScaLAPACK routines
 * wrappers for other numeric types can be added in the future
 */
template <typename number>
inline void
Cgerv2d(int /*context*/,
        int /*M*/,
        int /*N*/,
        number * /*A*/,
        int /*lda*/,
        int /*rsrc*/,
        int /*csrc*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
Cgerv2d(int context, int M, int N, double *A, int lda, int rsrc, int csrc)
{
  Cdgerv2d(context, M, N, A, lda, rsrc, csrc);
}

inline void
Cgerv2d(int context, int M, int N, float *A, int lda, int rsrc, int csrc)
{
  Csgerv2d(context, M, N, A, lda, rsrc, csrc);
}


template <typename number>
inline void
Cgesd2d(int /*context*/,
        int /*M*/,
        int /*N*/,
        number * /*A*/,
        int /*lda*/,
        int /*rdest*/,
        int /*cdest*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
Cgesd2d(int context, int M, int N, double *A, int lda, int rdest, int cdest)
{
  Cdgesd2d(context, M, N, A, lda, rdest, cdest);
}

inline void
Cgesd2d(int context, int M, int N, float *A, int lda, int rdest, int cdest)
{
  Csgesd2d(context, M, N, A, lda, rdest, cdest);
}


template <typename number>
inline void
ppotrf(const char * /*UPLO*/,
       const int * /*N*/,
       number * /*A*/,
       const int * /*IA*/,
       const int * /*JA*/,
       const int * /*DESCA*/,
       int * /*INFO*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
ppotrf(const char *UPLO,
       const int  *N,
       double     *A,
       const int  *IA,
       const int  *JA,
       const int  *DESCA,
       int        *INFO)
{
  pdpotrf_(UPLO, N, A, IA, JA, DESCA, INFO);
}

inline void
ppotrf(const char *UPLO,
       const int  *N,
       float      *A,
       const int  *IA,
       const int  *JA,
       const int  *DESCA,
       int        *INFO)
{
  pspotrf_(UPLO, N, A, IA, JA, DESCA, INFO);
}


template <typename number>
inline void
pgetrf(const int * /*m*/,
       const int * /*n*/,
       number * /*A*/,
       const int * /*IA*/,
       const int * /*JA*/,
       const int * /*DESCA*/,
       int * /*ipiv*/,
       int * /*INFO*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
pgetrf(const int *m,
       const int *n,
       double    *A,
       const int *IA,
       const int *JA,
       const int *DESCA,
       int       *ipiv,
       int       *INFO)
{
  pdgetrf_(m, n, A, IA, JA, DESCA, ipiv, INFO);
}

inline void
pgetrf(const int *m,
       const int *n,
       float     *A,
       const int *IA,
       const int *JA,
       const int *DESCA,
       int       *ipiv,
       int       *INFO)
{
  psgetrf_(m, n, A, IA, JA, DESCA, ipiv, INFO);
}


template <typename number>
inline void
ppotri(const char * /*UPLO*/,
       const int * /*N*/,
       number * /*A*/,
       const int * /*IA*/,
       const int * /*JA*/,
       const int * /*DESCA*/,
       int * /*INFO*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
ppotri(const char *UPLO,
       const int  *N,
       double     *A,
       const int  *IA,
       const int  *JA,
       const int  *DESCA,
       int        *INFO)
{
  pdpotri_(UPLO, N, A, IA, JA, DESCA, INFO);
}

inline void
ppotri(const char *UPLO,
       const int  *N,
       float      *A,
       const int  *IA,
       const int  *JA,
       const int  *DESCA,
       int        *INFO)
{
  pspotri_(UPLO, N, A, IA, JA, DESCA, INFO);
}


template <typename number>
inline void
pgetri(const int * /*N*/,
       number * /*A*/,
       const int * /*IA*/,
       const int * /*JA*/,
       const int * /*DESCA*/,
       const int * /*ipiv*/,
       number * /*work*/,
       int * /*lwork*/,
       int * /*iwork*/,
       int * /*liwork*/,
       int * /*info*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
pgetri(const int *N,
       double    *A,
       const int *IA,
       const int *JA,
       const int *DESCA,
       const int *ipiv,
       double    *work,
       int       *lwork,
       int       *iwork,
       int       *liwork,
       int       *info)
{
  pdgetri_(N, A, IA, JA, DESCA, ipiv, work, lwork, iwork, liwork, info);
}

inline void
pgetri(const int *N,
       float     *A,
       const int *IA,
       const int *JA,
       const int *DESCA,
       const int *ipiv,
       float     *work,
       int       *lwork,
       int       *iwork,
       int       *liwork,
       int       *info)
{
  psgetri_(N, A, IA, JA, DESCA, ipiv, work, lwork, iwork, liwork, info);
}

template <typename number>
inline void
ptrtri(const char * /*UPLO*/,
       const char * /*DIAG*/,
       const int * /*N*/,
       number * /*A*/,
       const int * /*IA*/,
       const int * /*JA*/,
       const int * /*DESCA*/,
       int * /*INFO*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
ptrtri(const char *UPLO,
       const char *DIAG,
       const int  *N,
       double     *A,
       const int  *IA,
       const int  *JA,
       const int  *DESCA,
       int        *INFO)
{
  pdtrtri_(UPLO, DIAG, N, A, IA, JA, DESCA, INFO);
}

inline void
ptrtri(const char *UPLO,
       const char *DIAG,
       const int  *N,
       float      *A,
       const int  *IA,
       const int  *JA,
       const int  *DESCA,
       int        *INFO)
{
  pstrtri_(UPLO, DIAG, N, A, IA, JA, DESCA, INFO);
}

template <typename number>
inline void
ppocon(const char * /*uplo*/,
       const int * /*N*/,
       const number * /*A*/,
       const int * /*IA*/,
       const int * /*JA*/,
       const int * /*DESCA*/,
       const number * /*ANORM*/,
       number * /*RCOND*/,
       number * /*WORK*/,
       const int * /*LWORK*/,
       int * /*IWORK*/,
       const int * /*LIWORK*/,
       int * /*INFO*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
ppocon(const char   *uplo,
       const int    *N,
       const double *A,
       const int    *IA,
       const int    *JA,
       const int    *DESCA,
       const double *ANORM,
       double       *RCOND,
       double       *WORK,
       const int    *LWORK,
       int          *IWORK,
       const int    *LIWORK,
       int          *INFO)
{
  pdpocon_(
    uplo, N, A, IA, JA, DESCA, ANORM, RCOND, WORK, LWORK, IWORK, LIWORK, INFO);
}

inline void
ppocon(const char  *uplo,
       const int   *N,
       const float *A,
       const int   *IA,
       const int   *JA,
       const int   *DESCA,
       const float *ANORM,
       float       *RCOND,
       float       *WORK,
       const int   *LWORK,
       int         *IWORK,
       const int   *LIWORK,
       int         *INFO)
{
  pspocon_(
    uplo, N, A, IA, JA, DESCA, ANORM, RCOND, WORK, LWORK, IWORK, LIWORK, INFO);
}


template <typename number>
inline number
plansy(const char * /*norm*/,
       const char * /*uplo*/,
       const int * /*N*/,
       const number * /*A*/,
       const int * /*IA*/,
       const int * /*JA*/,
       const int * /*DESCA*/,
       number * /*work*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline double
plansy(const char   *norm,
       const char   *uplo,
       const int    *N,
       const double *A,
       const int    *IA,
       const int    *JA,
       const int    *DESCA,
       double       *work)
{
  return pdlansy_(norm, uplo, N, A, IA, JA, DESCA, work);
}

inline float
plansy(const char  *norm,
       const char  *uplo,
       const int   *N,
       const float *A,
       const int   *IA,
       const int   *JA,
       const int   *DESCA,
       float       *work)
{
  return pslansy_(norm, uplo, N, A, IA, JA, DESCA, work);
}


template <typename number>
inline void
pgesv(const int * /*n*/,
      const int * /*nrhs*/,
      number * /*A*/,
      const int * /*ia*/,
      const int * /*ja*/,
      const int * /*desca*/,
      int * /*ipiv*/,
      number * /*B*/,
      const int * /*ib*/,
      const int * /*jb*/,
      const int * /*descb*/,
      int * /*info*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
pgesv(const int *n,
      const int *nrhs,
      double    *A,
      const int *ia,
      const int *ja,
      const int *desca,
      int       *ipiv,
      double    *B,
      const int *ib,
      const int *jb,
      const int *descb,
      int       *info)
{
  pdgesv_(n, nrhs, A, ia, ja, desca, ipiv, B, ib, jb, descb, info);
}

inline void
pgesv(const int *n,
      const int *nrhs,
      float     *A,
      const int *ia,
      const int *ja,
      const int *desca,
      int       *ipiv,
      float     *B,
      const int *ib,
      const int *jb,
      const int *descb,
      int       *info)
{
  psgesv_(n, nrhs, A, ia, ja, desca, ipiv, B, ib, jb, descb, info);
}


template <typename number>
inline void
pgemm(const char * /*transa*/,
      const char * /*transb*/,
      const int * /*m*/,
      const int * /*n*/,
      const int * /*k*/,
      const number * /*alpha*/,
      number * /*A*/,
      const int * /*IA*/,
      const int * /*JA*/,
      const int * /*DESCA*/,
      number * /*B*/,
      const int * /*IB*/,
      const int * /*JB*/,
      const int * /*DESCB*/,
      const number * /*beta*/,
      number * /*C*/,
      const int * /*IC*/,
      const int * /*JC*/,
      const int * /*DESCC*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
pgemm(const char   *transa,
      const char   *transb,
      const int    *m,
      const int    *n,
      const int    *k,
      const double *alpha,
      const double *A,
      const int    *IA,
      const int    *JA,
      const int    *DESCA,
      const double *B,
      const int    *IB,
      const int    *JB,
      const int    *DESCB,
      const double *beta,
      double       *C,
      const int    *IC,
      const int    *JC,
      const int    *DESCC)
{
  pdgemm_(transa,
          transb,
          m,
          n,
          k,
          alpha,
          A,
          IA,
          JA,
          DESCA,
          B,
          IB,
          JB,
          DESCB,
          beta,
          C,
          IC,
          JC,
          DESCC);
}

inline void
pgemm(const char  *transa,
      const char  *transb,
      const int   *m,
      const int   *n,
      const int   *k,
      const float *alpha,
      const float *A,
      const int   *IA,
      const int   *JA,
      const int   *DESCA,
      const float *B,
      const int   *IB,
      const int   *JB,
      const int   *DESCB,
      const float *beta,
      float       *C,
      const int   *IC,
      const int   *JC,
      const int   *DESCC)
{
  psgemm_(transa,
          transb,
          m,
          n,
          k,
          alpha,
          A,
          IA,
          JA,
          DESCA,
          B,
          IB,
          JB,
          DESCB,
          beta,
          C,
          IC,
          JC,
          DESCC);
}


template <typename number>
inline number
plange(const char * /*norm*/,
       const int * /*m*/,
       const int * /*n*/,
       const number * /*A*/,
       const int * /*ia*/,
       const int * /*ja*/,
       const int * /*desca*/,
       number * /*work*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline double
plange(const char   *norm,
       const int    *m,
       const int    *n,
       const double *A,
       const int    *ia,
       const int    *ja,
       const int    *desca,
       double       *work)
{
  return pdlange_(norm, m, n, A, ia, ja, desca, work);
}

inline float
plange(const char  *norm,
       const int   *m,
       const int   *n,
       const float *A,
       const int   *ia,
       const int   *ja,
       const int   *desca,
       float       *work)
{
  return pslange_(norm, m, n, A, ia, ja, desca, work);
}


template <typename number>
inline void
psyev(const char * /*jobz*/,
      const char * /*uplo*/,
      const int * /*m*/,
      number * /*A*/,
      const int * /*ia*/,
      const int * /*ja*/,
      int * /*desca*/,
      number * /*w*/,
      number * /*z*/,
      const int * /*iz*/,
      const int * /*jz*/,
      int * /*descz*/,
      number * /*work*/,
      const int * /*lwork*/,
      int * /*info*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
psyev(const char *jobz,
      const char *uplo,
      const int  *m,
      double     *A,
      const int  *ia,
      const int  *ja,
      int        *desca,
      double     *w,
      double     *z,
      const int  *iz,
      const int  *jz,
      int        *descz,
      double     *work,
      const int  *lwork,
      int        *info)
{
  pdsyev_(
    jobz, uplo, m, A, ia, ja, desca, w, z, iz, jz, descz, work, lwork, info);
}

inline void
psyev(const char *jobz,
      const char *uplo,
      const int  *m,
      float      *A,
      const int  *ia,
      const int  *ja,
      int        *desca,
      float      *w,
      float      *z,
      const int  *iz,
      const int  *jz,
      int        *descz,
      float      *work,
      const int  *lwork,
      int        *info)
{
  pssyev_(
    jobz, uplo, m, A, ia, ja, desca, w, z, iz, jz, descz, work, lwork, info);
}


template <typename number>
inline void
placpy(const char * /*uplo*/,
       const int * /*m*/,
       const int * /*n*/,
       const number * /*A*/,
       const int * /*ia*/,
       const int * /*ja*/,
       const int * /*desca*/,
       number * /*B*/,
       const int * /*ib*/,
       const int * /*jb*/,
       const int * /*descb*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
placpy(const char   *uplo,
       const int    *m,
       const int    *n,
       const double *A,
       const int    *ia,
       const int    *ja,
       const int    *desca,
       double       *B,
       const int    *ib,
       const int    *jb,
       const int    *descb)
{
  pdlacpy_(uplo, m, n, A, ia, ja, desca, B, ib, jb, descb);
}

inline void
placpy(const char  *uplo,
       const int   *m,
       const int   *n,
       const float *A,
       const int   *ia,
       const int   *ja,
       const int   *desca,
       float       *B,
       const int   *ib,
       const int   *jb,
       const int   *descb)
{
  pslacpy_(uplo, m, n, A, ia, ja, desca, B, ib, jb, descb);
}


template <typename number>
inline void
pgemr2d(const int * /*m*/,
        const int * /*n*/,
        const number * /*A*/,
        const int * /*ia*/,
        const int * /*ja*/,
        const int * /*desca*/,
        number * /*B*/,
        const int * /*ib*/,
        const int * /*jb*/,
        const int * /*descb*/,
        const int * /*ictxt*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
pgemr2d(const int    *m,
        const int    *n,
        const double *A,
        const int    *ia,
        const int    *ja,
        const int    *desca,
        double       *B,
        const int    *ib,
        const int    *jb,
        const int    *descb,
        const int    *ictxt)
{
  pdgemr2d_(m, n, A, ia, ja, desca, B, ib, jb, descb, ictxt);
}

inline void
pgemr2d(const int   *m,
        const int   *n,
        const float *A,
        const int   *ia,
        const int   *ja,
        const int   *desca,
        float       *B,
        const int   *ib,
        const int   *jb,
        const int   *descb,
        const int   *ictxt)
{
  psgemr2d_(m, n, A, ia, ja, desca, B, ib, jb, descb, ictxt);
}


template <typename number>
inline void
plamch(const int * /*ictxt*/, const char * /*cmach*/, number & /*val*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
plamch(const int *ictxt, const char *cmach, double &val)
{
  val = pdlamch_(ictxt, cmach);
}

inline void
plamch(const int *ictxt, const char *cmach, float &val)
{
  val = pslamch_(ictxt, cmach);
}


template <typename number>
inline void
psyevx(const char * /*jobz*/,
       const char * /*range*/,
       const char * /*uplo*/,
       const int * /*n*/,
       number * /*A*/,
       const int * /*ia*/,
       const int * /*ja*/,
       const int * /*desca*/,
       number * /*VL*/,
       number * /*VU*/,
       const int * /*il*/,
       const int * /*iu*/,
       number * /*abstol*/,
       const int * /*m*/,
       const int * /*nz*/,
       number * /*w*/,
       number * /*orfac*/,
       number * /*Z*/,
       const int * /*iz*/,
       const int * /*jz*/,
       const int * /*descz*/,
       number * /*work*/,
       int * /*lwork*/,
       int * /*iwork*/,
       int * /*liwork*/,
       int * /*ifail*/,
       int * /*iclustr*/,
       number * /*gap*/,
       int * /*info*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
psyevx(const char *jobz,
       const char *range,
       const char *uplo,
       const int  *n,
       double     *A,
       const int  *ia,
       const int  *ja,
       const int  *desca,
       double     *VL,
       double     *VU,
       const int  *il,
       const int  *iu,
       double     *abstol,
       const int  *m,
       const int  *nz,
       double     *w,
       double     *orfac,
       double     *Z,
       const int  *iz,
       const int  *jz,
       const int  *descz,
       double     *work,
       int        *lwork,
       int        *iwork,
       int        *liwork,
       int        *ifail,
       int        *iclustr,
       double     *gap,
       int        *info)
{
  pdsyevx_(jobz,
           range,
           uplo,
           n,
           A,
           ia,
           ja,
           desca,
           VL,
           VU,
           il,
           iu,
           abstol,
           m,
           nz,
           w,
           orfac,
           Z,
           iz,
           jz,
           descz,
           work,
           lwork,
           iwork,
           liwork,
           ifail,
           iclustr,
           gap,
           info);
}

inline void
psyevx(const char *jobz,
       const char *range,
       const char *uplo,
       const int  *n,
       float      *A,
       const int  *ia,
       const int  *ja,
       const int  *desca,
       float      *VL,
       float      *VU,
       const int  *il,
       const int  *iu,
       float      *abstol,
       const int  *m,
       const int  *nz,
       float      *w,
       float      *orfac,
       float      *Z,
       const int  *iz,
       const int  *jz,
       const int  *descz,
       float      *work,
       int        *lwork,
       int        *iwork,
       int        *liwork,
       int        *ifail,
       int        *iclustr,
       float      *gap,
       int        *info)
{
  pssyevx_(jobz,
           range,
           uplo,
           n,
           A,
           ia,
           ja,
           desca,
           VL,
           VU,
           il,
           iu,
           abstol,
           m,
           nz,
           w,
           orfac,
           Z,
           iz,
           jz,
           descz,
           work,
           lwork,
           iwork,
           liwork,
           ifail,
           iclustr,
           gap,
           info);
}


template <typename number>
inline void
pgesvd(const char * /*jobu*/,
       const char * /*jobvt*/,
       const int * /*m*/,
       const int * /*n*/,
       number * /*A*/,
       const int * /*ia*/,
       const int * /*ja*/,
       const int * /*desca*/,
       number * /*S*/,
       number * /*U*/,
       const int * /*iu*/,
       const int * /*ju*/,
       const int * /*descu*/,
       number * /*VT*/,
       const int * /*ivt*/,
       const int * /*jvt*/,
       const int * /*descvt*/,
       number * /*work*/,
       int * /*lwork*/,
       int * /*info*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
pgesvd(const char *jobu,
       const char *jobvt,
       const int  *m,
       const int  *n,
       double     *A,
       const int  *ia,
       const int  *ja,
       const int  *desca,
       double     *S,
       double     *U,
       const int  *iu,
       const int  *ju,
       const int  *descu,
       double     *VT,
       const int  *ivt,
       const int  *jvt,
       const int  *descvt,
       double     *work,
       int        *lwork,
       int        *info)
{
  pdgesvd_(jobu,
           jobvt,
           m,
           n,
           A,
           ia,
           ja,
           desca,
           S,
           U,
           iu,
           ju,
           descu,
           VT,
           ivt,
           jvt,
           descvt,
           work,
           lwork,
           info);
}

inline void
pgesvd(const char *jobu,
       const char *jobvt,
       const int  *m,
       const int  *n,
       float      *A,
       const int  *ia,
       const int  *ja,
       const int  *desca,
       float      *S,
       float      *U,
       const int  *iu,
       const int  *ju,
       const int  *descu,
       float      *VT,
       const int  *ivt,
       const int  *jvt,
       const int  *descvt,
       float      *work,
       int        *lwork,
       int        *info)
{
  psgesvd_(jobu,
           jobvt,
           m,
           n,
           A,
           ia,
           ja,
           desca,
           S,
           U,
           iu,
           ju,
           descu,
           VT,
           ivt,
           jvt,
           descvt,
           work,
           lwork,
           info);
}


template <typename number>
inline void
pgels(const char * /*trans*/,
      const int * /*m*/,
      const int * /*n*/,
      const int * /*nrhs*/,
      number * /*A*/,
      const int * /*ia*/,
      const int * /*ja*/,
      const int * /*desca*/,
      number * /*B*/,
      const int * /*ib*/,
      const int * /*jb*/,
      const int * /*descb*/,
      number * /*work*/,
      int * /*lwork*/,
      int * /*info*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
pgels(const char *trans,
      const int  *m,
      const int  *n,
      const int  *nrhs,
      double     *A,
      const int  *ia,
      const int  *ja,
      const int  *desca,
      double     *B,
      const int  *ib,
      const int  *jb,
      const int  *descb,
      double     *work,
      int        *lwork,
      int        *info)
{
  pdgels_(
    trans, m, n, nrhs, A, ia, ja, desca, B, ib, jb, descb, work, lwork, info);
}

inline void
pgels(const char *trans,
      const int  *m,
      const int  *n,
      const int  *nrhs,
      float      *A,
      const int  *ia,
      const int  *ja,
      const int  *desca,
      float      *B,
      const int  *ib,
      const int  *jb,
      const int  *descb,
      float      *work,
      int        *lwork,
      int        *info)
{
  psgels_(
    trans, m, n, nrhs, A, ia, ja, desca, B, ib, jb, descb, work, lwork, info);
}


template <typename number>
inline void
pgeadd(const char * /*transa*/,
       const int * /*m*/,
       const int * /*n*/,
       const number * /*alpha*/,
       const number * /*A*/,
       const int * /*IA*/,
       const int * /*JA*/,
       const int * /*DESCA*/,
       const number * /*beta*/,
       number * /*C*/,
       const int * /*IC*/,
       const int * /*JC*/,
       const int * /*DESCC*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
pgeadd(const char   *transa,
       const int    *m,
       const int    *n,
       const double *alpha,
       const double *A,
       const int    *IA,
       const int    *JA,
       const int    *DESCA,
       const double *beta,
       double       *C,
       const int    *IC,
       const int    *JC,
       const int    *DESCC)
{
  pdgeadd_(transa, m, n, alpha, A, IA, JA, DESCA, beta, C, IC, JC, DESCC);
}

inline void
pgeadd(const char  *transa,
       const int   *m,
       const int   *n,
       const float *alpha,
       const float *A,
       const int   *IA,
       const int   *JA,
       const int   *DESCA,
       const float *beta,
       float       *C,
       const int   *IC,
       const int   *JC,
       const int   *DESCC)
{
  psgeadd_(transa, m, n, alpha, A, IA, JA, DESCA, beta, C, IC, JC, DESCC);
}


template <typename number>
inline void
ptran(const int * /*m*/,
      const int * /*n*/,
      const number * /*alpha*/,
      const number * /*A*/,
      const int * /*IA*/,
      const int * /*JA*/,
      const int * /*DESCA*/,
      const number * /*beta*/,
      number * /*C*/,
      const int * /*IC*/,
      const int * /*JC*/,
      const int * /*DESCC*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
ptran(const int    *m,
      const int    *n,
      const double *alpha,
      const double *A,
      const int    *IA,
      const int    *JA,
      const int    *DESCA,
      const double *beta,
      double       *C,
      const int    *IC,
      const int    *JC,
      const int    *DESCC)
{
  pdtran_(m, n, alpha, A, IA, JA, DESCA, beta, C, IC, JC, DESCC);
}

inline void
ptran(const int   *m,
      const int   *n,
      const float *alpha,
      const float *A,
      const int   *IA,
      const int   *JA,
      const int   *DESCA,
      const float *beta,
      float       *C,
      const int   *IC,
      const int   *JC,
      const int   *DESCC)
{
  pstran_(m, n, alpha, A, IA, JA, DESCA, beta, C, IC, JC, DESCC);
}


template <typename number>
inline void
psyevr(const char * /*jobz*/,
       const char * /*range*/,
       const char * /*uplo*/,
       const int * /*n*/,
       number * /*A*/,
       const int * /*IA*/,
       const int * /*JA*/,
       const int * /*DESCA*/,
       const number * /*VL*/,
       const number * /*VU*/,
       const int * /*IL*/,
       const int * /*IU*/,
       int * /*m*/,
       int * /*nz*/,
       number * /*w*/,
       number * /*Z*/,
       const int * /*IZ*/,
       const int * /*JZ*/,
       const int * /*DESCZ*/,
       number * /*work*/,
       int * /*lwork*/,
       int * /*iwork*/,
       int * /*liwork*/,
       int * /*info*/)
{
  DEAL_II_NOT_IMPLEMENTED();
}

inline void
psyevr(const char   *jobz,
       const char   *range,
       const char   *uplo,
       const int    *n,
       double       *A,
       const int    *IA,
       const int    *JA,
       const int    *DESCA,
       const double *VL,
       const double *VU,
       const int    *IL,
       const int    *IU,
       int          *m,
       int          *nz,
       double       *w,
       double       *Z,
       const int    *IZ,
       const int    *JZ,
       const int    *DESCZ,
       double       *work,
       int          *lwork,
       int          *iwork,
       int          *liwork,
       int          *info)
{
  /*
   * Netlib ScaLAPACK performs floating point tests (e.g. divide-by-zero) within
   * the call to pdsyevr causing floating point exceptions to be thrown (at
   * least in debug mode). Therefore, we wrap the calls to pdsyevr into the
   * following code to suppress the exception.
   */
#  ifdef DEAL_II_HAVE_FP_EXCEPTIONS
  fenv_t fp_exceptions;
  feholdexcept(&fp_exceptions);
#  endif

  pdsyevr_(jobz,
           range,
           uplo,
           n,
           A,
           IA,
           JA,
           DESCA,
           VL,
           VU,
           IL,
           IU,
           m,
           nz,
           w,
           Z,
           IZ,
           JZ,
           DESCZ,
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
psyevr(const char  *jobz,
       const char  *range,
       const char  *uplo,
       const int   *n,
       float       *A,
       const int   *IA,
       const int   *JA,
       const int   *DESCA,
       const float *VL,
       const float *VU,
       const int   *IL,
       const int   *IU,
       int         *m,
       int         *nz,
       float       *w,
       float       *Z,
       const int   *IZ,
       const int   *JZ,
       const int   *DESCZ,
       float       *work,
       int         *lwork,
       int         *iwork,
       int         *liwork,
       int         *info)
{
  /*
   * Netlib ScaLAPACK performs floating point tests (e.g. divide-by-zero) within
   * the call to pssyevr causing floating point exceptions to be thrown (at
   * least in debug mode). Therefore, we wrap the calls to pssyevr into the
   * following code to suppress the exception.
   */
#  ifdef DEAL_II_HAVE_FP_EXCEPTIONS
  fenv_t fp_exceptions;
  feholdexcept(&fp_exceptions);
#  endif

  pssyevr_(jobz,
           range,
           uplo,
           n,
           A,
           IA,
           JA,
           DESCA,
           VL,
           VU,
           IL,
           IU,
           m,
           nz,
           w,
           Z,
           IZ,
           JZ,
           DESCZ,
           work,
           lwork,
           iwork,
           liwork,
           info);

#  ifdef DEAL_II_HAVE_FP_EXCEPTIONS
  fesetenv(&fp_exceptions);
#  endif
}

#endif // DEAL_II_WITH_SCALAPACK

// This file does not actually import anything into namespace dealii,
// but to avoid it being completely empty to some of our scripts, we
// need to make sure it opens and closes the namespace at least once.
DEAL_II_NAMESPACE_OPEN // Do not convert for module purposes
  DEAL_II_NAMESPACE_CLOSE


#endif // dealii_scalapack_templates_h
