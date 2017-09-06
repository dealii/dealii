// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

#include <deal.II/lac/scalapack.h>
#include <deal.II/base/mpi.h>

// useful examples:
// https://stackoverflow.com/questions/14147705/cholesky-decomposition-scalapack-error/14203864
// http://icl.cs.utk.edu/lapack-forum/viewtopic.php?t=139   // second post by Julien Langou
// https://andyspiros.wordpress.com/2011/07/08/an-example-of-blacs-with-c/
// http://qboxcode.org/trac/browser/qb/tags/rel1_63_4/src/Matrix.C
// https://gitlab.phys.ethz.ch/lwossnig/lecture/blob/a534f562dfb2ad5c564abe5c2356d5d956fb7218/examples/mpi/scalapack.cpp
// https://github.com/elemental/Elemental/blob/master/src/core/imports/scalapack.cpp
// https://scicomp.stackexchange.com/questions/7766/performance-optimization-or-tuning-possible-for-scalapack-gemm
//
// info:
// http://www.netlib.org/scalapack/slug/index.html       // User guide
// http://www.netlib.org/scalapack/slug/node135.html // How to Measure Errors

// FIXME: similar to lapack_templates.h move those to scalapack_template.h
extern "C"
{
  /* Basic Linear Algebra Communication Subprograms (BLACS) declarations */
  // https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dinitb.htm#dinitb
  void Cblacs_pinfo(int *, int *);
  void Cblacs_get(int icontxt, int what, int *val);
  void Cblacs_gridinit(int *context, const char *order, int grid_height, int grid_width);
  void Cblacs_gridinfo(int  context, int *grid_height, int *grid_width, int *grid_row, int *grid_col);
  void Cblacs_pcoord(int, int, int *, int *);
  void Cblacs_gridexit(int context);
  void Cblacs_barrier(int, const char *);
  void Cblacs_exit(int error_code);
  void Cdgerv2d(int, int, int, double *, int, int, int);
  void Cdgesd2d(int, int, int, double *, int, int, int);

  int Csys2blacs_handle(MPI_Comm comm);

  /**
   * NUMber of Rows Or Columns) -- computes how many rows and columns each process owns.
   *
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dnumy.htm
   */
  int numroc_ (const int *n, const int *nb, const int *iproc, const int *isproc, const int *nprocs);

  /**
   * Computes the Cholesky factorization of an N-by-N real
   *  symmetric positive definite distributed matrix sub( A ) denoting
   *  A(IA:IA+N-1, JA:JA+N-1).
   * see http://www.netlib.org/scalapack/explore-html/d5/d9e/pdpotrf_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lpotrf.htm
   */
  void pdpotrf_(const char *UPLO,
                const int *N,
                double *A, const int *IA, const int *JA, const int *DESCA,
                int *INFO);

  /**
   *  Computes the inverse of a real symmetric positive definite
   *  distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1) using the
   *  Cholesky factorization sub( A ) = U**T*U or L*L**T computed by
   *  PDPOTRF.
   *
   * see http://www.netlib.org/scalapack/explore-html/d2/d44/pdpotri_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lpotri.htm
   * https://software.intel.com/en-us/mkl-developer-reference-c-p-potri
   */
  void pdpotri_(const char *UPLO,
                const int *N,
                double *A, const int *IA, const int *JA, const int *DESCA,
                int *INFO);

  /**
   *  Estimates the reciprocal of the condition number (in the
   *  1-norm) of a real symmetric positive definite distributed matrix
   *  using the Cholesky factorization.
   *
   *  https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lpocon.htm#lpocon
   *  http://www.netlib.org/scalapack/explore-html/d4/df7/pdpocon_8f.html
   *  https://software.intel.com/en-us/mkl-developer-reference-fortran-pocon
   */
  void pdpocon_(const char *uplo,
                const int *N,
                const double *A, const int *IA, const int *JA, const int *DESCA,
                const double *ANORM, double *RCOND,
                double *WORK, const int *LWORK,
                int *IWORK, const int *LIWORK,
                int *INFO);

  /**
   * Norm of a real symmetric matrix
   *
   * http://www.netlib.org/scalapack/explore-html/dd/d12/pdlansy_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_pdlansy.htm#pdlansy
   */
  double pdlansy_(const char *norm,
                  const char *uplo,
                  const int *N,
                  const double *A, const int *IA, const int *JA, const int *DESCA,
                  double *work);

  /**
   *  Computes the Least Common Multiple (LCM) of two positive integers @p M and @p N.
   *  In fact the routine computes the greatest common divisor (GCD) and
   *  use the fact that M*N = GCD*LCM.
   *
   *  http://www.netlib.org/scalapack/explore-html/d0/d9b/ilcm_8f_source.html
   */
  int ilcm_(const int *M, const int *N);

  /**
   * returns the ceiling of the division of two integers.
   *
   * http://www.netlib.org/scalapack/explore-html/df/d07/iceil_8f_source.html
   */
  int iceil_(const int *i1, const int *i2);

  void descinit_ (int *desc, const int *m, const int *n, const int *mb, const int *nb, const int *irsrc, const int *icsrc, const int *ictxt, const int *lld, int *info);

  /**
   * computes the global index of a distributed matrix entry
   * pointed to by the local index @p indxloc of the process indicated by
   * @p iproc.
   *
   * @param indxloc The local index of the distributed matrix entry.
   * @param nb Block size, size of the blocks the distributed matrix is split into.
   * @param iproc The coordinate of the process whose local array row or column is to be determined
   * @param isrcproc  The coordinate of the process that possesses the first row/column of the distributed matrix
   * @param nprocs The total number processes over which the distributed matrix is distributed
   */
  int indxl2g_ (const int *indxloc, const int *nb, const int *iproc, const int *isrcproc, const int *nprocs);

  void pdgesv_(const int *n, const int *nrhs,
               double *A, const int *ia, const int *ja, const int *desca,
               int *ipiv,
               double *B, const int *ib, const int *jb, const int *descb,
               int *info);

  void pdgemm_(const char *transa, const char *transb,
               const int *m, const int *n, const int *k,
               const double *alpha,
               double *A, const int *IA, const int *JA, const int *DESCA,
               double *B, const int *IB, const int *JB, const int *DESCB,
               const double *beta,
               double *C, const int *IC, const int *JC, const int *DESCC);

  double pdlange_(char const *norm,
                  int const &m, int const &n,
                  double *A, int const &ia, int const &ja, int *desca,
                  double *work);

  int indxg2p_(const int *glob, const int *nb, const int *iproc, const int *isproc, const int *nprocs);

}

DEAL_II_NAMESPACE_OPEN

namespace
{
  /**
   * Internal function to create processor grid based on the total
   * number of cores and MPI size. See
   *
   * Amesos heuristics:
   * https://github.com/trilinos/Trilinos/blob/master/packages/amesos/src/Amesos_Scalapack.cpp#L142-L166
   *
   * Elemental default grid: El::Grid::Grid(mpi_comm,...)
   * https://github.com/elemental/Elemental/blob/master/src/core/Grid.cpp#L67-L91
   */
  inline
  void choose_the_processor_grid(int &n_process_rows, int &n_process_columns, const unsigned int size, const unsigned int block_size, const int n_processes)
  {
    // Few notes from the ScaLAPACK user guide:
    // It is possible to predict the best grid shape given the number of processes available:
    // Pr x Pc <= P
    // This, however, depends on the task to be done.
    // LU , QR and QL factorizations perform better for “flat” process grids (Pr < Pc )
    // For large N, Pc = 2*Pr is a good choice, whereas for small N, one should choose small Pr
    // Square or near square grids are more optimal for Cholesky factorization.
    // LQ and RQ factorizations take advantage of “tall” grids (Pr > Pc )

    // Below we always try to create 2D processor grids:

    // Get the total number of cores we can occupy in a square dense matrix
    // with square blocks when every core owns only a single block:
    const int n_processes_heuristic = int(std::ceil((1.*size)/block_size))*
                                      int(std::ceil((1.*size)/block_size));
    const int Np = std::min(n_processes_heuristic, n_processes);

    // Now we need to split Np into  Pr x Pc. Assume we know the shape/ratio
    // Pc =: ratio * Pr
    // therefore
    // Np = Pc * Pc / ratio
    const double ratio = 1.0;
    int Pc = std::sqrt(ratio * Np);

    // one could rounds up Pc to the number which has zero remainder from the division of Np
    // while ( Np % Pc != 0 )
    //  ++Pc;
    // but this affects the grid shape dramatically, i.e. 10 cores 3x3 becomes 2x5.

    // limit our estimate to be in [2, Np]
    n_process_columns = std::min (Np, std::max(2, Pc));
    // finally, get the rows:
    n_process_rows = Np / n_process_columns ;

    Assert (n_process_columns >=1 && n_process_rows >=1 && n_processes >= n_process_rows*n_process_columns,
            ExcMessage("error in process grid: "+
                       std::to_string(n_process_rows)+"x"+
                       std::to_string(n_process_columns)+
                       "="+
                       std::to_string(n_process_rows*n_process_columns)+
                       " out of " +
                       std::to_string(n_processes)));

    // For example,
    // 320x320 with 32x32 blocks and 16 cores:
    // Pc = 1.0 * Pr  => 4x4 grid
    // Pc = 0.5 * Pr  => 8x2 grid
    // Pc = 2.0 * Pr  => 3x5 grid
  }
}



template <typename NumberType>
ScaLAPACKMatrix<NumberType>::ScaLAPACKMatrix(const size_type size,
                                             MPI_Comm mpi_communicator,
                                             const unsigned int block_size)
  :
  TransposeTable<NumberType> (),
  state (LAPACKSupport::unusable),
  properties(LAPACKSupport::general),
  mpi_communicator(mpi_communicator),
  n_rows(size),
  n_columns(size),
  row_block_size(block_size),
  column_block_size(block_size),
  uplo('L'),
  first_process_row(0),
  first_process_column(0),
  submatrix_row(1),
  submatrix_column(1)
{
  // NOTE: routines that exploit symmetry (Cholesky pdpotrf), require row_block_size = column_block_size
  Assert (block_size > 0,
          ExcMessage("Block size has to be positive."));
  Assert (block_size <= size,
          ExcMessage("BLock size can not be greater than the size of the matrix"));

  MPI_Comm_size(mpi_communicator, &n_mpi_processes);
  MPI_Comm_rank(mpi_communicator, &this_mpi_process);

  choose_the_processor_grid(n_process_rows, n_process_columns, size, block_size, n_mpi_processes);

  // processor grid order
  // FIXME: default to column major
  const bool column_major = false;

  // Initialize Cblas context from the provided communicator
  blacs_context = Csys2blacs_handle(mpi_communicator);
  const char *order = ( column_major ? "Col" : "Row" );
  // FIXME: blacs_context can be modified below. Thus Cblacs2sys_handle
  // may not return the same MPI communicator
  Cblacs_gridinit(&blacs_context, order, n_process_rows, n_process_columns);

  // Blacs may modify the grid size on processes which are not used
  // in the grid. So provide copies below:
  int procrows_ = n_process_rows;
  int proccols_ = n_process_columns;
  Cblacs_gridinfo( blacs_context, &procrows_, &proccols_, &this_process_row, &this_process_column );

  // If this MPI core is not on the grid, flag is as inactive and
  // skip all jobs
  // FIXME: different condition is used here
  // https://stackoverflow.com/questions/18516915/calling-blacs-with-more-processes-than-used
  if (this_process_row < 0 || this_process_column < 0)
    active = false;
  else
    active = true;

  if (active)
    {
      Assert (n_process_rows == proccols_ && proccols_ == n_process_columns,
              ExcInternalError());

      // Get local sizes:
      n_local_rows = numroc_(&n_rows, &row_block_size, &this_process_row, &first_process_row, &n_process_rows);
      n_local_columns = numroc_(&n_columns, &column_block_size, &this_process_column, &first_process_column, &n_process_columns);

      // LLD_A = MAX(1,NUMROC(M_A, MB_A, MYROW, RSRC_A, NPROW)), different between processes
      int lda = std::max(1,n_local_rows);

      int info=0;
      descinit_(descriptor, &n_rows, &n_columns, &row_block_size, &column_block_size,&first_process_row,&first_process_column,&blacs_context, &lda, &info);
      AssertThrow (info==0, LAPACKSupport::ExcErrorCode("descinit_", info));

      this->reinit(n_local_rows, n_local_columns);
    }
  else
    {
      // set process-local variables to something telling:
      n_local_rows = -1;
      n_local_columns = -1;
      for (unsigned int i = 0; i < 9; ++i)
        descriptor[i] = -1;
    }

  // Create an auxiliary communicator which has root and all inactive cores
  // Assume that inactive cores start with id=n_process_rows*n_process_columns
  Assert (active || this_mpi_process >= n_process_rows*n_process_columns,
          ExcInternalError());

  std::vector<int> inactive_with_root_ranks;
  inactive_with_root_ranks.push_back(0);
  for (int i = n_process_rows*n_process_columns; i < n_mpi_processes; ++i)
    inactive_with_root_ranks.push_back(i);

  // Get the group of processes in mpi_communicator
  int ierr = 0;
  MPI_Group all_group;
  ierr = MPI_Comm_group(mpi_communicator, &all_group);
  AssertThrowMPI(ierr);

  // Construct the group containing all ranks we need:
  MPI_Group inactive_with_root_group;
  const int n = inactive_with_root_ranks.size();
  ierr = MPI_Group_incl(all_group,
                        n, &inactive_with_root_ranks[0],
                        &inactive_with_root_group);
  AssertThrowMPI(ierr);

  // Create the communicator based on the group
  // Note that, on most cores the communicator will be MPI_COMM_NULL
  // FIXME: switch to MPI_Comm_create_group for MPI-3 so that only processes within the to-be subgroup call this
  ierr = MPI_Comm_create(mpi_communicator, inactive_with_root_group,
                         &mpi_communicator_inactive_with_root);
  AssertThrowMPI(ierr);

  MPI_Group_free(&all_group);
  MPI_Group_free(&inactive_with_root_group);

  // Double check that the process with rank 0 in subgroup is active:
#ifdef DEBUG
  if (mpi_communicator_inactive_with_root != MPI_COMM_NULL)
    {
      int subgroup_rank = -1, subgroup_size = -1;
      MPI_Comm_rank(mpi_communicator_inactive_with_root, &subgroup_rank);
      MPI_Comm_size(mpi_communicator_inactive_with_root, &subgroup_size);
      if (subgroup_rank == 0)
        Assert (active, ExcInternalError());
    }
#endif
}



template <typename NumberType>
ScaLAPACKMatrix<NumberType> &
ScaLAPACKMatrix<NumberType>::operator = (const FullMatrix<NumberType> &matrix)
{
  // FIXME: another way to copy is to use pdgeadd_ PBLAS routine.
  // This routine computes sum of two matrices B:=a*A+b*B.
  // Matrices can have different distribution,in particular matrixA can
  // be owned by only one process, so we can set a=1 and b=0 to copy
  // non-distributed matrix A into distributed matrix B.
  Assert (n_rows == int(matrix.m()), ExcDimensionMismatch(n_rows, matrix.m()));
  Assert (n_columns == int(matrix.n()), ExcDimensionMismatch(n_columns, matrix.n()));

  if (active)
    {
      for (int i=0; i < n_local_rows; ++i)
        {
          const int glob_i = global_row(i);
          for (int j = 0; j < n_local_columns; ++j)
            {
              const int glob_j = global_column(j);
              local_el(i,j) = matrix(glob_i, glob_j);
            }
        }
    }
  // FIXME: move it from operator = to copy_from() with a const bool symmetric flag.
  properties = LAPACKSupport::symmetric;
  state = LAPACKSupport::matrix;
  return *this;
}



template <typename NumberType>
int
ScaLAPACKMatrix<NumberType>::m() const
{
  return n_rows;
}



template <typename NumberType>
int
ScaLAPACKMatrix<NumberType>::n() const
{
  return n_columns;
}



template <typename NumberType>
int
ScaLAPACKMatrix<NumberType>::local_m() const
{
  return n_local_rows;
}



template <typename NumberType>
int
ScaLAPACKMatrix<NumberType>::local_n() const
{
  return n_local_columns;
}



template <typename NumberType>
int
ScaLAPACKMatrix<NumberType>::global_row(const int loc_row) const
{
  Assert (loc_row >= 0 && loc_row < n_local_rows,
          ExcIndexRange(loc_row,0,n_local_rows));
  const int i = loc_row+1;
  return indxl2g_ (&i, &row_block_size, &this_process_row, &first_process_row, &n_process_rows) - 1;
}



template <typename NumberType>
int
ScaLAPACKMatrix<NumberType>::global_column(const int loc_column) const
{
  Assert (loc_column >= 0 && loc_column < n_local_columns,
          ExcIndexRange(loc_column,0,n_local_columns));
  const int j = loc_column+1;
  return indxl2g_ (&j, &column_block_size, &this_process_column, &first_process_column, &n_process_columns) - 1;
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::copy_to (FullMatrix<NumberType> &matrix) const
{
  // FIXME: use PDGEMR2D for copying?
  // PDGEMR2D copies a submatrix of A on a submatrix of B.
  // A and B can have different distributions
  // see http://icl.cs.utk.edu/lapack-forum/viewtopic.php?t=50
  Assert (n_rows == int(matrix.m()), ExcDimensionMismatch(n_rows, matrix.m()));
  Assert (n_columns == int(matrix.n()), ExcDimensionMismatch(n_columns, matrix.n()));

  if (active)
    {
      matrix = 0.;
      for (int i=0; i < n_local_rows; ++i)
        {
          const int glob_i = global_row(i);
          for (int j = 0; j < n_local_columns; ++j)
            {
              const int glob_j = global_column(j);
              matrix(glob_i, glob_j) = local_el(i,j);
            }
        }
    }

  Utilities::MPI::sum(matrix, mpi_communicator, matrix);

  // we could move the following lines under the main loop above,
  // but they would be dependent on glob_i and glob_j, which
  // won't make it much prettier
  if (properties == LAPACKSupport::lower_triangular)
    for (unsigned int i = 0; i < matrix.n(); ++i)
      for (unsigned int j = i+1; j < matrix.m(); ++j)
        matrix(i,j) = (state == LAPACKSupport::inverse_matrix ? matrix(j,i) : 0.);
  else if (properties == LAPACKSupport::upper_triangular)
    for (unsigned int i = 0; i < matrix.n(); ++i)
      for (unsigned int j = 0; j < i; ++j)
        matrix(i,j) = (state == LAPACKSupport::inverse_matrix ? matrix(j,i) : 0.);
}



template <typename NumberType>
ScaLAPACKMatrix<NumberType>::~ScaLAPACKMatrix()
{
  if (active)
    Cblacs_gridexit(blacs_context);

  MPI_Comm_free(&mpi_communicator_inactive_with_root);
}



template <typename NumberType>
void ScaLAPACKMatrix<NumberType>::compute_cholesky_factorization()
{
  Assert (n_columns == n_rows,
          ExcMessage("Cholesky factorization can be applied to SPD matrices only."));
  if (active)
    {
      int info = 0;
      NumberType *A_loc = &this->values[0];
      pdpotrf_(&uplo,&n_columns,A_loc,&submatrix_row,&submatrix_column,descriptor,&info);
      AssertThrow (info==0, LAPACKSupport::ExcErrorCode("pdpotrf", info));
    }
  properties = (uplo=='L' ? LAPACKSupport::lower_triangular : LAPACKSupport::upper_triangular);
  state = LAPACKSupport::cholesky;
}



template <typename NumberType>
void ScaLAPACKMatrix<NumberType>::invert()
{
  if (state == LAPACKSupport::matrix)
    compute_cholesky_factorization();

  if (active)
    {
      int info = 0;
      NumberType *A_loc = &this->values[0];
      pdpotri_ (&uplo,&n_columns, A_loc, &submatrix_row, &submatrix_column, descriptor,&info);
      AssertThrow (info==0, LAPACKSupport::ExcErrorCode("pdpotri", info));
    }
  state = LAPACKSupport::inverse_matrix;
}



template <typename NumberType>
NumberType ScaLAPACKMatrix<NumberType>::reciprocal_condition_number(const NumberType a_norm) const
{
  Assert (state == LAPACKSupport::cholesky,
          ExcMessage("Matrix has to be in Cholesky state before calling this function."));
  NumberType rcond = 0.;
  if (active)
    {
      int lwork = 2 * n_local_rows + 3 * n_local_columns + column_block_size;
      int liwork = n_local_rows;
      work.resize(lwork);
      iwork.resize(liwork);
      int info = 0;
      const NumberType *A_loc = &this->values[0];
      pdpocon_(&uplo, &n_columns, A_loc, &submatrix_row, &submatrix_column, descriptor,
               &a_norm, &rcond, &work[0], &lwork, &iwork[0], &liwork, &info);
      AssertThrow (info==0, LAPACKSupport::ExcErrorCode("pdpocon", info));
    }
  send_to_inactive(rcond);
  return rcond;
}



template <typename NumberType>
NumberType ScaLAPACKMatrix<NumberType>::l1_norm() const
{
  const char type('O');
  return norm(type);
}



template <typename NumberType>
NumberType ScaLAPACKMatrix<NumberType>::linfty_norm() const
{
  const char type('I');
  return norm(type);
}



template <typename NumberType>
NumberType ScaLAPACKMatrix<NumberType>::frobenius_norm() const
{
  const char type('F');
  return norm(type);
}



template <typename NumberType>
NumberType ScaLAPACKMatrix<NumberType>::norm(const char type) const
{
  Assert (state == LAPACKSupport::matrix ||
          state == LAPACKSupport::inverse_matrix,
          ExcMessage("norms can be called in matrix state only."));

  NumberType res = 0.;

  if (active)
    {
      //int IROFFA = MOD( IA-1, MB_A )
      //int ICOFFA = MOD( JA-1, NB_A )
      const int lcm = ilcm_(&n_process_rows, &n_process_columns);
      const int v2 = lcm/n_process_rows;

      const int IAROW = indxg2p_(&submatrix_row, &row_block_size, &this_process_row, &first_process_row, &n_process_rows);
      const int IACOL = indxg2p_(&submatrix_column, &column_block_size, &this_process_column, &first_process_column, &n_process_columns);
      const int Np0   = numroc_(&n_columns/*+IROFFA*/, &row_block_size, &this_process_row, &IAROW, &n_process_rows);
      const int Nq0   = numroc_(&n_columns/*+ICOFFA*/, &column_block_size, &this_process_column, &IACOL, &n_process_columns);

      const int v1 = iceil_(&Np0, &row_block_size);
      const int ldw = (n_local_rows==n_local_columns) ?
                      0 :
                      row_block_size*iceil_(&v1,&v2);

      const int lwork = (type == 'M' || type == 'F' || type == 'E' ) ?
                        0 :
                        2*Nq0+Np0+ldw;
      work.resize(lwork);
      const NumberType *A_loc = &this->values[0];
      res = pdlansy_(&type, &uplo, &n_columns, A_loc, &submatrix_row, &submatrix_column, descriptor, &work[0]);
    }
  send_to_inactive(res);
  return res;
}



template <typename NumberType>
void ScaLAPACKMatrix<NumberType>::send_to_inactive(NumberType &value) const
{
  if (mpi_communicator_inactive_with_root != MPI_COMM_NULL)
    {
      MPI_Bcast(&value,1,MPI_DOUBLE,
                0/*from root*/,
                mpi_communicator_inactive_with_root);
    }
}



template <typename NumberType>
int ScaLAPACKMatrix<NumberType>::get_process_grid_rows() const
{
  return n_process_rows;
}



template <typename NumberType>
int ScaLAPACKMatrix<NumberType>::get_process_grid_columns() const
{
  return n_process_columns;
}



// instantiations
template class ScaLAPACKMatrix<double>;

DEAL_II_NAMESPACE_CLOSE
