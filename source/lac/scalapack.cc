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
                                             const std::shared_ptr<const Utilities::MPI::ProcessGrid> &process_grid,
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
ScaLAPACKMatrix<NumberType>::global_row(const unsigned int loc_row) const
{
  Assert (n_local_rows >= 0 && loc_row < static_cast<unsigned int>(n_local_rows),
          ExcIndexRange(loc_row,0,n_local_rows));
  const int i = loc_row+1;
  return indxl2g_ (&i, &row_block_size, &(grid->this_process_row), &first_process_row, &(grid->n_process_rows)) - 1;
}



template <typename NumberType>
unsigned int
ScaLAPACKMatrix<NumberType>::global_column(const unsigned int loc_column) const
{
  Assert (n_local_columns >= 0 && loc_column < static_cast<unsigned int>(n_local_columns),
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
void
ScaLAPACKMatrix<NumberType>::copy_to (ScaLAPACKMatrix<NumberType> &dest) const
{
  Assert (n_rows == dest.n_rows, ExcDimensionMismatch(n_rows, dest.n_rows));
  Assert (n_columns == dest.n_columns, ExcDimensionMismatch(n_columns, dest.n_columns));

  if (this->grid->mpi_process_is_active)
    AssertThrow (this->descriptor[0]==1,ExcMessage("Copying of ScaLAPACK matrices only implemented for dense matrices"));
  if (dest.grid->mpi_process_is_active)
    AssertThrow (dest.descriptor[0]==1,ExcMessage("Copying of ScaLAPACK matrices only implemented for dense matrices"));

  /*
   * just in case of different process grids or block-cyclic distributions
   * inter-process communication is necessary
   * if distributed matrices have the same process grid and block sizes, local copying is enough
   */
  if ( (this->grid != dest.grid) || (row_block_size != dest.row_block_size) || (column_block_size != dest.column_block_size) )
    {
      /*
       * get the MPI communicator, which is the union of the source and destination MPI communicator
       */
      int ierr = 0;
      MPI_Group group_source, group_dest, group_union;
      ierr = MPI_Comm_group(this->grid->mpi_communicator, &group_source);
      AssertThrowMPI(ierr);
      ierr = MPI_Comm_group(dest.grid->mpi_communicator, &group_dest);
      AssertThrowMPI(ierr);
      ierr = MPI_Group_union(group_source, group_dest, &group_union);
      AssertThrowMPI(ierr);
      MPI_Comm mpi_communicator_union;
      // to create a communicator representing the union of the source and destination MPI communicator
      // we need a communicator containing all  desired processes --> use MPI_COMM_WORLD
      ierr = MPI_Comm_create_group(MPI_COMM_WORLD, group_union, 5, &mpi_communicator_union);
      AssertThrowMPI(ierr);

      /*
       * The routine pgemr2d requires a BLACS context resembling at least the union of process grids
       * described by the BLACS contexts of the source and destination matrix
       */
      int union_blacs_context = Csys2blacs_handle(mpi_communicator_union);
      const char *order = "Col";
      int union_n_process_rows = Utilities::MPI::n_mpi_processes(mpi_communicator_union);
      int union_n_process_columns = 1;
      Cblacs_gridinit(&union_blacs_context, order, union_n_process_rows, union_n_process_columns);

      const NumberType *loc_vals_source = NULL;
      NumberType *loc_vals_dest = NULL;

      if (this->grid->mpi_process_is_active && (this->values.size()>0))
        {
          AssertThrow(this->values.size()>0,dealii::ExcMessage("source: process is active but local matrix empty"));
          loc_vals_source = &this->values[0];
        }
      if (dest.grid->mpi_process_is_active && (dest.values.size()>0))
        {
          AssertThrow(dest.values.size()>0,dealii::ExcMessage("destination: process is active but local matrix empty"));
          loc_vals_dest = &dest.values[0];
        }
      pgemr2d(&n_rows, &n_columns, loc_vals_source, &submatrix_row, &submatrix_column, descriptor,
              loc_vals_dest, &dest.submatrix_row, &dest.submatrix_column, dest.descriptor,
              &union_blacs_context);

      Cblacs_gridexit(union_blacs_context);

      if (mpi_communicator_union != MPI_COMM_NULL)
        MPI_Comm_free(&mpi_communicator_union);
      MPI_Group_free(&group_source);
      MPI_Group_free(&group_dest);
      MPI_Group_free(&group_union);
    }
  else
    //process is active in the process grid
    if (this->grid->mpi_process_is_active)
      dest.values = this->values;

  dest.state = state;
  dest.property = property;
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
      //pdpotrf_(&uplo,&n_columns,A_loc,&submatrix_row,&submatrix_column,descriptor,&info);
      ppotrf(&uplo,&n_columns,A_loc,&submatrix_row,&submatrix_column,descriptor,&info);
      AssertThrow (info==0, LAPACKSupport::ExcErrorCode("ppotrf", info));
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
      ppotri (&uplo,&n_columns, A_loc, &submatrix_row, &submatrix_column, descriptor,&info);
      AssertThrow (info==0, LAPACKSupport::ExcErrorCode("ppotri", info));
    }
  state = LAPACKSupport::inverse_matrix;
}



template <typename NumberType>
std::vector<NumberType> ScaLAPACKMatrix<NumberType>::eigenpairs_symmetric_by_index(const std::pair<unsigned int,unsigned int> &index_limits,
    const bool compute_eigenvectors)
{
  // check validity of index limits
  Assert (index_limits.first < (unsigned int)n_rows,ExcIndexRange(index_limits.first,0,n_rows));
  Assert (index_limits.second < (unsigned int)n_rows,ExcIndexRange(index_limits.second,0,n_rows));

  std::pair<unsigned int,unsigned int> idx = std::make_pair(std::min(index_limits.first,index_limits.second),
                                                            std::max(index_limits.first,index_limits.second));

  // compute all eigenvalues/eigenvectors
  if (idx.first==0 && idx.second==(unsigned int)n_rows-1)
    return eigenpairs_symmetric(compute_eigenvectors);
  else
    return eigenpairs_symmetric(compute_eigenvectors,idx);
}



template <typename NumberType>
std::vector<NumberType> ScaLAPACKMatrix<NumberType>::eigenpairs_symmetric_by_value(const std::pair<NumberType,NumberType> &value_limits,
    const bool compute_eigenvectors)
{
  Assert (!std::isnan(value_limits.first),ExcMessage("value_limits.first is NaN"));
  Assert (!std::isnan(value_limits.second),ExcMessage("value_limits.second is NaN"));

  std::pair<unsigned int,unsigned int> indices = std::make_pair(numbers::invalid_unsigned_int,numbers::invalid_unsigned_int);

  return eigenpairs_symmetric(compute_eigenvectors,indices,value_limits);
}



template <typename NumberType>
std::vector<NumberType>
ScaLAPACKMatrix<NumberType>::eigenpairs_symmetric(const bool compute_eigenvectors,
                                                  const std::pair<unsigned int, unsigned int> &eigenvalue_idx,
                                                  const std::pair<NumberType,NumberType> &eigenvalue_limits)
{
  Assert (state == LAPACKSupport::matrix,
          ExcMessage("Matrix has to be in Matrix state before calling this function."));
  Assert (property == LAPACKSupport::symmetric,
          ExcMessage("Matrix has to be symmetric for this operation."));

  Threads::Mutex::ScopedLock lock (mutex);

  const bool use_values = (std::isnan(eigenvalue_limits.first) || std::isnan(eigenvalue_limits.second)) ? false : true;
  const bool use_indices = ((eigenvalue_idx.first==numbers::invalid_unsigned_int) || (eigenvalue_idx.second==numbers::invalid_unsigned_int)) ? false : true;

  Assert(!(use_values && use_indices),ExcMessage("Prescribing both the index and value range for the eigenvalues is ambiguous"));

  // if computation of eigenvectors is not required use a sufficiently small distributed matrix
  std::unique_ptr<ScaLAPACKMatrix<NumberType>> eigenvectors = compute_eigenvectors ?
                                                              std::make_unique<ScaLAPACKMatrix<NumberType>>(n_rows,grid,row_block_size) :
                                                              std::make_unique<ScaLAPACKMatrix<NumberType>>(grid->n_process_rows,grid->n_process_columns,grid,1,1);

  eigenvectors->property = property;
  // number of eigenvalues to be returned; upon successful exit ev contains the m seclected eigenvalues in ascending order
  int m = n_rows;
  std::vector<NumberType> ev(n_rows);

  if (grid->mpi_process_is_active)
    {
      int info = 0;
      /*
       * for jobz==N only eigenvalues are computed, for jobz='V' also the eigenvectors of the matrix are computed
       */
      char jobz = compute_eigenvectors ? 'V' : 'N';
      char range;
      // default value is to compute all eigenvalues and optionally eigenvectors
      bool all_eigenpairs=true;
      NumberType vl,vu;
      int il,iu;
      // number of eigenvectors to be returned;
      // upon successful exit the first m=nz columns contain the selected eigenvectors (only if jobz=='V')
      int nz;
      NumberType abstol;
      char cmach = compute_eigenvectors ? 'U' : 'S';

      // orfac decides which eigenvectors should be reorthogonalized
      // see http://www.netlib.org/scalapack/explore-html/df/d1a/pdsyevx_8f_source.html for explanation
      // to keeps simple no reorthogonalized will be done by setting orfac to 0
      NumberType orfac = 0;
      //contains the indices of eigenvectors that failed to converge
      std::vector<int> ifail;
      // This array contains indices of eigenvectors corresponding to
      // a cluster of eigenvalues that could not be reorthogonalized
      // due to insufficient workspace
      // see http://www.netlib.org/scalapack/explore-html/df/d1a/pdsyevx_8f_source.html for explanation
      std::vector<int> iclustr;
      // This array contains the gap between eigenvalues whose
      // eigenvectors could not be reorthogonalized.
      // see http://www.netlib.org/scalapack/explore-html/df/d1a/pdsyevx_8f_source.html for explanation
      std::vector<NumberType> gap(n_local_rows * n_local_columns);

      // index range for eigenvalues is not specified
      if (!use_indices)
        {
          // interval for eigenvalues is not specified and consequently all eigenvalues/eigenpairs will be computed
          if (!use_values)
            {
              range = 'A';
              all_eigenpairs = true;
            }
          else
            {
              range = 'V';
              all_eigenpairs = false;
              vl = std::min(eigenvalue_limits.first,eigenvalue_limits.second);
              vu = std::max(eigenvalue_limits.first,eigenvalue_limits.second);
            }
        }
      else
        {
          range = 'I';
          all_eigenpairs = false;
          //as Fortran starts counting/indexing from 1 unlike C/C++, where it starts from 0
          il = std::min(eigenvalue_idx.first,eigenvalue_idx.second) + 1;
          iu = std::max(eigenvalue_idx.first,eigenvalue_idx.second) + 1;
        }
      NumberType *A_loc = &this->values[0];
      /*
       * by setting lwork to -1 a workspace query for optimal length of work is performed
       */
      int lwork=-1;
      int liwork=-1;
      NumberType *eigenvectors_loc = (compute_eigenvectors ? &eigenvectors->values[0] : nullptr);
      work.resize(1);
      iwork.resize (1);

      if (all_eigenpairs)
        {
          psyev(&jobz, &uplo, &n_rows, A_loc, &submatrix_row, &submatrix_column, descriptor, &ev[0],
                eigenvectors_loc, &eigenvectors->submatrix_row, &eigenvectors->submatrix_column, eigenvectors->descriptor,
                &work[0], &lwork, &info);
        }
      else
        {
          plamch( &(this->grid->blacs_context), &cmach, abstol);
          abstol *= 2;
          ifail.resize(n_rows);
          iclustr.resize(n_local_rows * n_local_columns);
          gap.resize(n_local_rows * n_local_columns);

          psyevx(&jobz, &range, &uplo, &n_rows, A_loc, &submatrix_row, &submatrix_column, descriptor,
                 &vl, &vu, &il, &iu, &abstol, &m, &nz, &ev[0], &orfac,
                 eigenvectors_loc, &eigenvectors->submatrix_row, &eigenvectors->submatrix_column, eigenvectors->descriptor,
                 &work[0], &lwork, &iwork[0], &liwork, &ifail[0], &iclustr[0], &gap[0], &info);
        }
      lwork=work[0];
      work.resize (lwork);

      if (all_eigenpairs)
        {
          psyev(&jobz, &uplo, &n_rows, A_loc, &submatrix_row, &submatrix_column, descriptor, &ev[0],
                eigenvectors_loc, &eigenvectors->submatrix_row, &eigenvectors->submatrix_column, eigenvectors->descriptor,
                &work[0], &lwork, &info);

          AssertThrow (info==0, LAPACKSupport::ExcErrorCode("psyev", info));
        }
      else
        {
          liwork = iwork[0];
          AssertThrow(liwork>0,ExcInternalError());
          iwork.resize(liwork);

          psyevx(&jobz, &range, &uplo, &n_rows, A_loc, &submatrix_row, &submatrix_column, descriptor,
                 &vl, &vu, &il, &iu, &abstol, &m, &nz, &ev[0], &orfac,
                 eigenvectors_loc, &eigenvectors->submatrix_row, &eigenvectors->submatrix_column, eigenvectors->descriptor,
                 &work[0], &lwork, &iwork[0], &liwork, &ifail[0], &iclustr[0], &gap[0], &info);

          AssertThrow (info==0, LAPACKSupport::ExcErrorCode("psyevx", info));
        }
      // if eigenvectors are queried copy eigenvectors to original matrix
      // as the temporary matrix eigenvectors has identical dimensions and
      // block-cyclic distribution we simply swap the local array
      if (compute_eigenvectors)
        this->values.swap(eigenvectors->values);

      //adapt the size of ev to fit m upon return
      while ((int)ev.size() > m)
        ev.pop_back();
    }
  /*
   * send number of computed eigenvalues to inactive processes
   */
  grid->send_to_inactive(&m, 1);

  /*
   * inactive processes have to resize array of eigenvalues
   */
  if (! grid->mpi_process_is_active)
    ev.resize (m);
  /*
   * send the eigenvalues to processors not being part of the process grid
   */
  grid->send_to_inactive(ev.data(), ev.size());

  /*
   * if only eigenvalues are queried the content of the matrix will be destroyed
   * if the eigenpairs are queried matrix A on exit stores the eigenvectors in the columns
   */
  if (compute_eigenvectors)
    {
      property = LAPACKSupport::Property::general;
      state = LAPACKSupport::eigenvalues;
    }
  else
    state = LAPACKSupport::unusable;

  return ev;
}



template <typename NumberType>
std::vector<NumberType> ScaLAPACKMatrix<NumberType>::compute_SVD(ScaLAPACKMatrix<NumberType> *U,
    ScaLAPACKMatrix<NumberType> *VT)
{
  Assert (state == LAPACKSupport::matrix,
          ExcMessage("Matrix has to be in Matrix state before calling this function."));
  Assert(row_block_size==column_block_size,ExcDimensionMismatch(row_block_size,column_block_size));

  const bool left_singluar_vectors = (U != nullptr) ? true : false;
  const bool right_singluar_vectors = (VT != nullptr) ? true : false;

  if (left_singluar_vectors)
    {
      Assert(n_rows==U->n_rows,ExcDimensionMismatch(n_rows,U->n_rows));
      Assert(U->n_rows==U->n_columns,ExcDimensionMismatch(U->n_rows,U->n_columns));
      Assert(row_block_size==U->row_block_size,ExcDimensionMismatch(row_block_size,U->row_block_size));
      Assert(column_block_size==U->column_block_size,ExcDimensionMismatch(column_block_size,U->column_block_size));
      Assert(grid->blacs_context==U->grid->blacs_context,ExcDimensionMismatch(grid->blacs_context,U->grid->blacs_context));
    }
  if (right_singluar_vectors)
    {
      Assert(n_columns==VT->n_rows,ExcDimensionMismatch(n_columns,VT->n_rows));
      Assert(VT->n_rows==VT->n_columns,ExcDimensionMismatch(VT->n_rows,VT->n_columns));
      Assert(row_block_size==VT->row_block_size,ExcDimensionMismatch(row_block_size,VT->row_block_size));
      Assert(column_block_size==VT->column_block_size,ExcDimensionMismatch(column_block_size,VT->column_block_size));
      Assert(grid->blacs_context==VT->grid->blacs_context,ExcDimensionMismatch(grid->blacs_context,VT->grid->blacs_context));
    }
  Threads::Mutex::ScopedLock lock (mutex);

  std::vector<NumberType> sv(std::min(n_rows,n_columns));

  if (grid->mpi_process_is_active)
    {
      char jobu = left_singluar_vectors ? 'V' : 'N';
      char jobvt = right_singluar_vectors ? 'V' : 'N';
      NumberType *A_loc = &this->values[0];
      NumberType *U_loc = left_singluar_vectors ? &(U->values[0]) : nullptr;
      NumberType *VT_loc = right_singluar_vectors ? &(VT->values[0]) : nullptr;
      int info = 0;
      /*
       * by setting lwork to -1 a workspace query for optimal length of work is performed
       */
      int lwork=-1;
      work.resize(1);

      pgesvd(&jobu,&jobvt,&n_rows,&n_columns,A_loc,&submatrix_row,&submatrix_column,descriptor,
             & *sv.begin(),U_loc,&U->submatrix_row,&U->submatrix_column,U->descriptor,
             VT_loc,&VT->submatrix_row,&VT->submatrix_column,VT->descriptor,
             &work[0],&lwork,&info);
      AssertThrow (info==0, LAPACKSupport::ExcErrorCode("pgesvd", info));

      lwork=work[0];
      work.resize(lwork);

      pgesvd(&jobu,&jobvt,&n_rows,&n_columns,A_loc,&submatrix_row,&submatrix_column,descriptor,
             & *sv.begin(),U_loc,&U->submatrix_row,&U->submatrix_column,U->descriptor,
             VT_loc,&VT->submatrix_row,&VT->submatrix_column,VT->descriptor,
             &work[0],&lwork,&info);
      AssertThrow (info==0, LAPACKSupport::ExcErrorCode("pgesvd", info));
    }

  /*
   * send the singular values to processors not being part of the process grid
   */
  grid->send_to_inactive(sv.data(), sv.size());

  property = LAPACKSupport::Property::general;
  state = LAPACKSupport::State::unusable;

  return sv;
}



template <typename NumberType>
void ScaLAPACKMatrix<NumberType>::least_squares(ScaLAPACKMatrix<NumberType> &B,
                                                const bool transpose)
{
  Assert(grid==B.grid,ExcMessage("The matrices A and B need to have the same process grid"));
  Assert (state == LAPACKSupport::matrix,
          ExcMessage("Matrix has to be in Matrix state before calling this function."));
  Assert (B.state == LAPACKSupport::matrix,
          ExcMessage("Matrix B has to be in Matrix state before calling this function."));

  transpose ?
  (Assert(n_columns==B.n_rows,ExcDimensionMismatch(n_columns,B.n_rows))) :
  (Assert(n_rows==B.n_rows,ExcDimensionMismatch(n_rows,B.n_rows)));

  //see https://www.ibm.com/support/knowledgecenter/en/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lgels.htm
  Assert(row_block_size==column_block_size,ExcMessage("Use identical block sizes for rows and columns of matrix A"));
  Assert(B.row_block_size==B.column_block_size,ExcMessage("Use identical block sizes for rows and columns of matrix B"));
  Assert(row_block_size==B.row_block_size,ExcMessage("Use identical block-cyclic distribution for matrices A and B"));

  Threads::Mutex::ScopedLock lock (mutex);

  if (grid->mpi_process_is_active)
    {
      char trans = transpose ? 'T' : 'N';
      NumberType *A_loc = & this->values[0];
      NumberType *B_loc = & B.values[0];
      int info = 0;
      /*
       * by setting lwork to -1 a workspace query for optimal length of work is performed
       */
      int lwork=-1;
      work.resize(1);

      pgels(&trans,&n_rows,&n_columns,&B.n_columns,A_loc,&submatrix_row,&submatrix_column,descriptor,
            B_loc,&B.submatrix_row,&B.submatrix_column,B.descriptor,&work[0],&lwork,&info);
      AssertThrow (info==0, LAPACKSupport::ExcErrorCode("pgels", info));

      lwork=work[0];
      work.resize(lwork);

      pgels(&trans,&n_rows,&n_columns,&B.n_columns,A_loc,&submatrix_row,&submatrix_column,descriptor,
            B_loc,&B.submatrix_row,&B.submatrix_column,B.descriptor,&work[0],&lwork,&info);
      AssertThrow (info==0, LAPACKSupport::ExcErrorCode("pgels", info));
    }
  state = LAPACKSupport::State::unusable;
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
      ppocon(&uplo, &n_columns, A_loc, &submatrix_row, &submatrix_column, descriptor,
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
      res = plansy(&type, &uplo, &n_columns, A_loc, &submatrix_row, &submatrix_column, descriptor, &work[0]);
    }
  grid->send_to_inactive(&res);
  return res;
}

// instantiations
template class ScaLAPACKMatrix<double>;
template class ScaLAPACKMatrix<float>;


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SCALAPACK
