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
#include <deal.II/lac/scalapack.templates.h>

DEAL_II_NAMESPACE_OPEN

template <typename NumberType>
ScaLAPACKMatrix<NumberType>::ScaLAPACKMatrix(const size_type n_rows_,
                                             const size_type n_columns_,
                                             const std::shared_ptr<const Utilities::MPI::ProcessGrid> process_grid,
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
                                             const std::shared_ptr<const Utilities::MPI::ProcessGrid> process_grid,
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


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SCALAPACK
