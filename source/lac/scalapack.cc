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

#ifdef DEAL_II_WITH_SCALAPACK

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>

#include <deal.II/base/conditional_ostream.h>

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
  void Cblacs_gridinit(int *context, const char *order, int grid_height, int grid_width);

  /**
   * Return the process row and column index.
   *
   * https://www.ibm.com/support/knowledgecenter/en/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dbinfo.htm
   */
  void Cblacs_gridinfo(int context, int *grid_height, int *grid_width, int *grid_row, int *grid_col);

  /**
   * Given the system process number, return the row and column coordinates in the BLACS' process grid.
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
  void Cdgerv2d(int context, int M, int N, double *A, int lda, int rsrc, int csrc);

  /**
   * Sends the general rectangular matrix A to the destination
   * process @p rdest @p cdest in the process grid.
   *
   * https://software.intel.com/en-us/mkl-developer-reference-c-2018-beta-gesd2d
   */
  void Cdgesd2d(int context , int M, int N, double *A, int lda, int rdest, int cdest);

  /**
   * Get BLACS context from MPI @p comm.
   */
  int Csys2blacs_handle(MPI_Comm comm);

  /**
   * Compute how many rows and columns each process owns (NUMber of Rows Or Columns).
   *
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dnumy.htm
   */
  int numroc_ (const int *n, const int *nb, const int *iproc, const int *isproc, const int *nprocs);

  /**
   * Compute the Cholesky factorization of an N-by-N real
   * symmetric positive definite distributed matrix sub( A ) denoting
   * A(IA:IA+N-1, JA:JA+N-1).
   *
   * http://www.netlib.org/scalapack/explore-html/d5/d9e/pdpotrf_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lpotrf.htm
   */
  void pdpotrf_(const char *UPLO,
                const int *N,
                double *A, const int *IA, const int *JA, const int *DESCA,
                int *INFO);

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
                const int *N,
                double *A, const int *IA, const int *JA, const int *DESCA,
                int *INFO);

  /**
   * Estimate the reciprocal of the condition number (in the
   * l1-norm) of a real symmetric positive definite distributed matrix
   * using the Cholesky factorization.
   *
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lpocon.htm#lpocon
   * http://www.netlib.org/scalapack/explore-html/d4/df7/pdpocon_8f.html
   * https://software.intel.com/en-us/mkl-developer-reference-fortran-pocon
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
  void descinit_ (int *desc,
                  const int *m, const int *n, const int *mb, const int *nb,
                  const int *irsrc, const int *icsrc,
                  const int *ictxt, const int *lld, int *info);

  /**
   * Compute the global index of a distributed matrix entry
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

  /**
   * Compute the solution to a real system of linear equations
   */
  void pdgesv_(const int *n, const int *nrhs,
               double *A, const int *ia, const int *ja, const int *desca,
               int *ipiv,
               double *B, const int *ib, const int *jb, const int *descb,
               int *info);

  /**
   * Perform one of the matrix-matrix operations:
   * sub( C ) := alpha*op( sub( A ) )*op( sub( B ) ) + beta*sub( C ),
   * where
   * sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),  and, op( X )  is one  of
   * op( X ) = X   or   op( X ) = X'.
   */
  void pdgemm_(const char *transa, const char *transb,
               const int *m, const int *n, const int *k,
               const double *alpha,
               double *A, const int *IA, const int *JA, const int *DESCA,
               double *B, const int *IB, const int *JB, const int *DESCB,
               const double *beta,
               double *C, const int *IC, const int *JC, const int *DESCC);

  /**
   * Return the value of the one norm, or the Frobenius norm, or the infinity norm,
   * or the element of largest absolute value of a distributed matrix
   */
  double pdlange_(char const *norm,
                  int const &m, int const &n,
                  double *A, int const &ia, int const &ja, int *desca,
                  double *work);

  /**
   * Compute the process coordinate which possesses the entry of a
   * distributed matrix specified by a global index
   */
  int indxg2p_(const int *glob, const int *nb, const int *iproc, const int *isproc, const int *nprocs);

  /**
   * Compute all eigenvalues and, optionally, eigenvectors of a real symmetric matrix A
   * by calling the recommended sequence of ScaLAPACK routines. In its present form, the routine assumes a homogeneous system
   * and makes no checks for consistency of the eigenvalues or eigenvectors across the different processes.
   * Because of this, it is possible that a heterogeneous system may return incorrect results without any error messages.
   *
   * http://www.netlib.org/scalapack/explore-html/d0/d1a/pdsyev_8f.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lsyev.htm#lsyev
   */
  void pdsyev_(const char *jobz, const char *uplo,
               const int *m, double *A, const int *ia, const int *ja, int *desca,
               double *w,
               double *z, const int *iz, const int *jz, int *descz,
               double *work, const int *lwork, int *info);

  /**
   * Copy all or a part of a distributed matrix A to another
   * distributed matrix B. No communication is performed, pdlacpy
   * performs a local copy sub(A) := sub(B), where sub(A) denotes
   * A(ia:ia+m-1,ja:ja+n-1) and sub(B) denotes B(ib:ib+m-1,jb:jb+n-1)
   */
  void pdlacpy_(const char *uplo,
                const int *m, const int *n, double *A, const int *ia, const int *ja, int *desca,
                double *B, const int *ib, const int *jb, int *descb);
}



DEAL_II_NAMESPACE_OPEN

namespace
{
  /**
   * Internal function to determine dimension of process grid based on the total
   * number of cores, matrix dimensions and matrix block sizes
   *
   * Amesos heuristics:
   * https://github.com/trilinos/Trilinos/blob/master/packages/amesos/src/Amesos_Scalapack.cpp#L142-L166
   *
   * Elemental default grid: El::Grid::Grid(mpi_comm,...)
   * https://github.com/elemental/Elemental/blob/master/src/core/Grid.cpp#L67-L91
   */
  inline
  std::pair<int,int> compute_processor_grid_sizes(MPI_Comm mpi_comm, const unsigned int m, const unsigned int n,
                                                  const unsigned int block_size_m, const unsigned int block_size_n)
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

    int n_processes;
    MPI_Comm_size(mpi_comm, &n_processes);

    // Get the total number of cores we can occupy in a rectangular dense matrix
    // with rectangular blocks when every core owns only a single block:
    const int n_processes_heuristic = int(std::ceil((1.*m)/block_size_m))*
                                      int(std::ceil((1.*n)/block_size_n));
    const int Np = std::min(n_processes_heuristic, n_processes);

    // Now we need to split Np into  Pr x Pc. Assume we know the shape/ratio
    // Pc =: ratio * Pr
    // therefore
    // Np = Pc * Pc / ratio
    // for quadratic matrices the ratio equals 1
    const double ratio = double(n)/m;
    int Pc = std::floor(std::sqrt(ratio * Np));

    // one could rounds up Pc to the number which has zero remainder from the division of Np
    // while ( Np % Pc != 0 )
    //  ++Pc;
    // but this affects the grid shape dramatically, i.e. 10 cores 3x3 becomes 2x5.

    // limit our estimate to be in [2, Np]
    int n_process_columns = std::min (Np, std::max(2, Pc));
    // finally, get the rows:
    int n_process_rows = Np / n_process_columns ;

    Assert (n_process_columns >=1 && n_process_rows >=1 && n_processes >= n_process_rows*n_process_columns,
            ExcMessage("error in process grid: "+
                       std::to_string(n_process_rows)+"x"+
                       std::to_string(n_process_columns)+
                       "="+
                       std::to_string(n_process_rows*n_process_columns)+
                       " out of " +
                       std::to_string(n_processes)));

    return std::make_pair(n_process_rows,n_process_columns);

    // For example,
    // 320x320 with 32x32 blocks and 16 cores:
    // Pc = 1.0 * Pr  => 4x4 grid
    // Pc = 0.5 * Pr  => 8x2 grid
    // Pc = 2.0 * Pr  => 3x5 grid
  }
}



ProcessGrid::ProcessGrid(MPI_Comm mpi_comm,
                         const std::pair<unsigned int,unsigned int> &grid_dimensions)
  :
  mpi_communicator(mpi_comm),
  this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator)),
  n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator)),
  n_process_rows(grid_dimensions.first),
  n_process_columns(grid_dimensions.second)
{
  Assert (grid_dimensions.first > 0,
          ExcMessage("Number of process grid rows has to be positive."));
  Assert (grid_dimensions.second > 0,
          ExcMessage("Number of process grid columns has to be positive."));

  Assert (grid_dimensions.first*grid_dimensions.second <= n_mpi_processes,
          ExcMessage("Size of process grid is larger than number of available MPI processes."));

  // processor grid order.
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

  // If this MPI core is not on the grid, flag it as inactive and
  // skip all jobs
  // FIXME: different condition is used here
  // https://stackoverflow.com/questions/18516915/calling-blacs-with-more-processes-than-used
  if (this_process_row < 0 || this_process_column < 0)
    mpi_process_is_active = false;
  else
    mpi_process_is_active = true;

  // Create an auxiliary communicator which has root and all inactive cores
  // Assume that inactive cores start with id=n_process_rows*n_process_columns
  Assert (mpi_process_is_active || this_mpi_process >= n_process_rows*n_process_columns,
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
                        n, inactive_with_root_ranks.data(),
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
  if (mpi_communicator_inactive_with_root != MPI_COMM_NULL &&
      Utilities::MPI::this_mpi_process(mpi_communicator_inactive_with_root) == 0)
    Assert (mpi_process_is_active, ExcInternalError());
#endif
}



ProcessGrid::ProcessGrid(MPI_Comm mpi_comm,
                         const unsigned int n_rows_matrix,
                         const unsigned int n_columns_matrix,
                         const unsigned int row_block_size,
                         const unsigned int column_block_size)
  :
  ProcessGrid(mpi_comm,
              compute_processor_grid_sizes(mpi_comm, n_rows_matrix, n_columns_matrix,
                                           row_block_size, column_block_size) )
{}



ProcessGrid::ProcessGrid(MPI_Comm mpi_comm,
                         const unsigned int n_rows,
                         const unsigned int n_columns)
  :
  ProcessGrid(mpi_comm,
              std::make_pair(n_rows,n_columns))
{}




ProcessGrid::~ProcessGrid()
{
  if (mpi_process_is_active)
    Cblacs_gridexit(blacs_context);

  MPI_Comm_free(&mpi_communicator_inactive_with_root);
}



unsigned int ProcessGrid::get_process_grid_rows() const
{
  return n_process_rows;
}



unsigned int ProcessGrid::get_process_grid_columns() const
{
  return n_process_columns;
}



bool ProcessGrid::is_process_active() const
{
  return mpi_process_is_active;
}



template <typename NumberType>
void ProcessGrid::send_to_inactive(NumberType *value, const int count) const
{
  Assert (count>0, ExcInternalError());
  if (mpi_communicator_inactive_with_root != MPI_COMM_NULL)
    {
      const int ierr =
        MPI_Bcast(value,count,
                  Utilities::MPI::internal::mpi_type_id (value),
                  0/*from root*/,
                  mpi_communicator_inactive_with_root);
      AssertThrowMPI(ierr);
    }
}



template <typename NumberType>
ScaLAPACKMatrix<NumberType>::ScaLAPACKMatrix(const size_type n_rows_,
                                             const size_type n_columns_,
                                             const std::shared_ptr<const ProcessGrid> process_grid,
                                             const size_type row_block_size_,
                                             const size_type column_block_size_,
                                             const LAPACKSupport::Property property)
  :
  TransposeTable<NumberType> (),
  state (LAPACKSupport::unusable),
  property(property),
  grid (process_grid),
  n_rows(n_rows_),
  n_columns(n_columns_),
  row_block_size(row_block_size_),
  column_block_size(column_block_size_),
  uplo('L'), // for non-symmetric matrices this is not needed
  first_process_row(0),
  first_process_column(0),
  submatrix_row(1),
  submatrix_column(1)
{
  Assert (row_block_size > 0,
          ExcMessage("Row block size has to be positive."));
  Assert (column_block_size > 0,
          ExcMessage("Column block size has to be positive."));
  Assert (row_block_size <= n_rows,
          ExcMessage("Row block size can not be greater than the number of rows of the matrix"));
  Assert (column_block_size <= n_columns,
          ExcMessage("Column block size can not be greater than the number of columns of the matrix"));

  if (grid->mpi_process_is_active)
    {
      // Get local sizes:
      n_local_rows = numroc_(&n_rows, &row_block_size, &(grid->this_process_row), &first_process_row, &(grid->n_process_rows));
      n_local_columns = numroc_(&n_columns, &column_block_size, &(grid->this_process_column), &first_process_column, &(grid->n_process_columns));

      // LLD_A = MAX(1,NUMROC(M_A, MB_A, MYROW, RSRC_A, NPROW)), different between processes
      int lda = std::max(1,n_local_rows);

      int info=0;
      descinit_(descriptor, &n_rows, &n_columns,
                &row_block_size, &column_block_size,
                &first_process_row, &first_process_column,
                &(grid->blacs_context), &lda, &info);
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
}



template <typename NumberType>
ScaLAPACKMatrix<NumberType>::ScaLAPACKMatrix(const size_type size,
                                             const std::shared_ptr<const ProcessGrid> process_grid,
                                             const size_type block_size,
                                             const LAPACKSupport::Property property)
  :
  ScaLAPACKMatrix<NumberType>(size,
                              size,
                              process_grid,
                              block_size,
                              block_size,
                              property)
{}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::set_property(const LAPACKSupport::Property property_)
{
  property = property_;
}



template <typename NumberType>
ScaLAPACKMatrix<NumberType> &
ScaLAPACKMatrix<NumberType>::operator = (const FullMatrix<NumberType> &matrix)
{
  // FIXME: another way to copy is to use pdgeadd_ PBLAS routine.
  // This routine computes the sum of two matrices B:=a*A+b*B.
  // Matrices can have different distribution,in particular matrix A can
  // be owned by only one process, so we can set a=1 and b=0 to copy
  // non-distributed matrix A into distributed matrix B.
  Assert (n_rows == int(matrix.m()), ExcDimensionMismatch(n_rows, matrix.m()));
  Assert (n_columns == int(matrix.n()), ExcDimensionMismatch(n_columns, matrix.n()));

  if (grid->mpi_process_is_active)
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
  state = LAPACKSupport::matrix;
  return *this;
}



template <typename NumberType>
unsigned int
ScaLAPACKMatrix<NumberType>::m() const
{
  return n_rows;
}



template <typename NumberType>
unsigned int
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
  return indxl2g_ (&i, &row_block_size, &(grid->this_process_row), &first_process_row, &(grid->n_process_rows)) - 1;
}



template <typename NumberType>
int
ScaLAPACKMatrix<NumberType>::global_column(const int loc_column) const
{
  Assert (loc_column >= 0 && loc_column < n_local_columns,
          ExcIndexRange(loc_column,0,n_local_columns));
  const int j = loc_column+1;
  return indxl2g_ (&j, &column_block_size, &(grid->this_process_column), &first_process_column, &(grid->n_process_columns)) - 1;
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

  if (grid->mpi_process_is_active)
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
  Utilities::MPI::sum(matrix, grid->mpi_communicator, matrix);

  // we could move the following lines under the main loop above,
  // but they would be dependent on glob_i and glob_j, which
  // won't make it much prettier
  if (property == LAPACKSupport::lower_triangular)
    for (unsigned int i = 0; i < matrix.n(); ++i)
      for (unsigned int j = i+1; j < matrix.m(); ++j)
        matrix(i,j) = (state == LAPACKSupport::inverse_matrix ? matrix(j,i) : 0.);
  else if (property == LAPACKSupport::upper_triangular)
    for (unsigned int i = 0; i < matrix.n(); ++i)
      for (unsigned int j = 0; j < i; ++j)
        matrix(i,j) = (state == LAPACKSupport::inverse_matrix ? matrix(j,i) : 0.);
}



template <typename NumberType>
void ScaLAPACKMatrix<NumberType>::compute_cholesky_factorization()
{
  Assert (n_columns == n_rows,
          ExcMessage("Cholesky factorization can be applied to SPD matrices only."));

  if (grid->mpi_process_is_active)
    {
      int info = 0;
      NumberType *A_loc = &this->values[0];
      pdpotrf_(&uplo,&n_columns,A_loc,&submatrix_row,&submatrix_column,descriptor,&info);
      AssertThrow (info==0, LAPACKSupport::ExcErrorCode("pdpotrf", info));
    }
  property = (uplo=='L' ? LAPACKSupport::lower_triangular : LAPACKSupport::upper_triangular);
  state = LAPACKSupport::cholesky;
}



template <typename NumberType>
void ScaLAPACKMatrix<NumberType>::invert()
{
  if (state == LAPACKSupport::matrix)
    compute_cholesky_factorization();

  if (grid->mpi_process_is_active)
    {
      int info = 0;
      NumberType *A_loc = &this->values[0];
      pdpotri_ (&uplo,&n_columns, A_loc, &submatrix_row, &submatrix_column, descriptor,&info);
      AssertThrow (info==0, LAPACKSupport::ExcErrorCode("pdpotri", info));
    }
  state = LAPACKSupport::inverse_matrix;
}



template <typename NumberType>
std::vector<NumberType> ScaLAPACKMatrix<NumberType>::eigenvalues_symmetric()
{
  Assert (state == LAPACKSupport::matrix,
          ExcMessage("Matrix has to be in Matrix state before calling this function."));
  Assert (property == LAPACKSupport::symmetric,
          ExcMessage("Matrix has to be symmetric for this operation."));
  Threads::Mutex::ScopedLock lock (mutex);

  ScaLAPACKMatrix<NumberType> Z (grid->n_mpi_processes, grid, 1);
  std::vector<NumberType> ev (n_rows);

  if (grid->mpi_process_is_active)
    {
      int info = 0;

      char jobz = 'N';
      NumberType *A_loc = &this->values[0];

      /*
       * by setting lwork to -1 a workspace query for optimal length of work is performed
      */
      int lwork=-1;
      NumberType *Z_loc = &Z.values[0];
      work.resize(1);

      pdsyev_(&jobz, &uplo, &n_rows, A_loc, &submatrix_row, &submatrix_column, descriptor, &ev[0],
              Z_loc, &Z.submatrix_row, &Z.submatrix_column, Z.descriptor, &work[0], &lwork, &info);

      lwork=work[0];
      work.resize (lwork);

      pdsyev_(&jobz, &uplo, &n_rows, A_loc, &submatrix_row, &submatrix_column, descriptor, &ev[0],
              Z_loc, &Z.submatrix_row, &Z.submatrix_column, Z.descriptor, &work[0], &lwork, &info);

      AssertThrow (info==0, LAPACKSupport::ExcErrorCode("pdsyev", info));
    }
  /*
   * send the eigenvalues to processors not being part of the process grid
   */
  grid->send_to_inactive(ev.data(), ev.size());

  /*
   * On exit, the lower triangle (if uplo='L') or the upper triangle (if uplo='U') of A,
   * including the diagonal, is destroyed. Therefore, the matrix is unusable
   */
  state = LAPACKSupport::unusable;

  return ev;
}



template <typename NumberType>
std::vector<NumberType> ScaLAPACKMatrix<NumberType>::eigenpairs_symmetric()
{
  Assert (state == LAPACKSupport::matrix,
          ExcMessage("Matrix has to be in Matrix state before calling this function."));
  Assert (property == LAPACKSupport::symmetric,
          ExcMessage("Matrix has to be symmetric for this operation."));

  Threads::Mutex::ScopedLock lock (mutex);

  ScaLAPACKMatrix<NumberType> eigenvectors (n_rows, grid, row_block_size);
  eigenvectors.property = property;
  std::vector<NumberType> ev(n_rows);

  if (grid->mpi_process_is_active)
    {
      int info = 0;

      /*
       * for jobz = 'V' all eigenpairs of the matrix are computed
       */
      char jobz = 'V';
      NumberType *A_loc = &this->values[0];

      /*
       * by setting lwork to -1 a workspace query for optimal length of work is performed
       */
      int lwork=-1;
      NumberType *eigenvectors_loc = &eigenvectors.values[0];
      work.resize(1);

      pdsyev_(&jobz, &uplo, &n_rows, A_loc, &submatrix_row, &submatrix_column, descriptor, &ev[0],
              eigenvectors_loc, &eigenvectors.submatrix_row, &eigenvectors.submatrix_column, eigenvectors.descriptor, &work[0], &lwork, &info);

      lwork=work[0];
      work.resize (lwork);

      pdsyev_(&jobz, &uplo, &n_rows, A_loc, &submatrix_row, &submatrix_column, descriptor, &ev[0],
              eigenvectors_loc, &eigenvectors.submatrix_row, &eigenvectors.submatrix_column, eigenvectors.descriptor, &work[0], &lwork, &info);

      AssertThrow (info==0, LAPACKSupport::ExcErrorCode("pdsyev", info));

      // copy eigenvectors to original matrix
      // as the temporary matrix eigenvectors has identical dimensions and
      // block-cyclic distribution we simply swap the local array
      this->values.swap(eigenvectors.values);
    }
  /*
   * send the eigenvalues to processors not being part of the process grid
   */
  grid->send_to_inactive(ev.data(), ev.size());

  /*
   *  On exit matrix A stores the eigenvectors in the columns
   */
  property = LAPACKSupport::Property::general;
  state = LAPACKSupport::eigenvalues;

  return ev;
}



template <typename NumberType>
NumberType ScaLAPACKMatrix<NumberType>::reciprocal_condition_number(const NumberType a_norm) const
{
  Assert (state == LAPACKSupport::cholesky,
          ExcMessage("Matrix has to be in Cholesky state before calling this function."));
  Threads::Mutex::ScopedLock lock (mutex);
  NumberType rcond = 0.;

  if (grid->mpi_process_is_active)
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
  grid->send_to_inactive(&rcond);
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
  Threads::Mutex::ScopedLock lock (mutex);
  NumberType res = 0.;

  if (grid->mpi_process_is_active)
    {
      //int IROFFA = MOD( IA-1, MB_A )
      //int ICOFFA = MOD( JA-1, NB_A )
      const int lcm = ilcm_(&(grid->n_process_rows), &(grid->n_process_columns));
      const int v2 = lcm/(grid->n_process_rows);

      const int IAROW = indxg2p_(&submatrix_row, &row_block_size, &(grid->this_process_row), &first_process_row, &(grid->n_process_rows));
      const int IACOL = indxg2p_(&submatrix_column, &column_block_size, &(grid->this_process_column), &first_process_column, &(grid->n_process_columns));
      const int Np0   = numroc_(&n_columns/*+IROFFA*/, &row_block_size, &(grid->this_process_row), &IAROW, &(grid->n_process_rows));
      const int Nq0   = numroc_(&n_columns/*+ICOFFA*/, &column_block_size, &(grid->this_process_column), &IACOL, &(grid->n_process_columns));

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
  grid->send_to_inactive(&res);
  return res;
}

// instantiations
template class ScaLAPACKMatrix<double>;
template void ProcessGrid::send_to_inactive<double>(double *, const int) const;


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SCALAPACK
