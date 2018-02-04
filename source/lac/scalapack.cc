// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

#include <deal.II/base/std_cxx14/memory.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>
#include <deal.II/lac/scalapack.templates.h>

#ifdef DEAL_II_WITH_HDF5
#include <hdf5.h>
#endif

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_HDF5

template<typename number>
inline hid_t hdf5_type_id (const number *)
{
  Assert (false, dealii::ExcNotImplemented());
  //don't know what to put here; it does not matter
  return -1;
}

inline hid_t hdf5_type_id (const double *)
{
  return H5T_NATIVE_DOUBLE;
}

inline hid_t hdf5_type_id (const float *)
{
  return H5T_NATIVE_FLOAT;
}

inline hid_t hdf5_type_id (const int *)
{
  return H5T_NATIVE_INT;
}

inline hid_t hdf5_type_id (const unsigned int *)
{
  return H5T_NATIVE_UINT;
}

inline hid_t hdf5_type_id (const char *)
{
  return H5T_NATIVE_CHAR;
}
#endif // DEAL_II_WITH_HDF5



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
ScaLAPACKMatrix<NumberType>::copy_to(ScaLAPACKMatrix<NumberType> &B,
                                     const std::pair<unsigned int,unsigned int> &offset_A,
                                     const std::pair<unsigned int,unsigned int> &offset_B,
                                     const std::pair<unsigned int,unsigned int> &submatrix_size) const
{
  //submatrix is empty
  if (submatrix_size.first == 0 || submatrix_size.second == 0)
    return;

  //range checking for matrix A
  Assert (offset_A.first<(unsigned int)(n_rows-submatrix_size.first+1),
          ExcIndexRange(offset_A.first,0,n_rows-submatrix_size.first+1));
  Assert (offset_A.second<(unsigned int)(n_columns-submatrix_size.second+1),
          ExcIndexRange(offset_A.second,0,n_columns-submatrix_size.second+1));

  //range checking for matrix B
  Assert (offset_B.first<(unsigned int)(B.n_rows-submatrix_size.first+1),
          ExcIndexRange(offset_B.first,0,B.n_rows-submatrix_size.first+1));
  Assert (offset_B.second<(unsigned int)(B.n_columns-submatrix_size.second+1),
          ExcIndexRange(offset_B.second,0,B.n_columns-submatrix_size.second+1));

  //Currently, copying of matrices will only be supported if A and B share the same MPI communicator
  int ierr, comparison;
  ierr = MPI_Comm_compare(grid->mpi_communicator,B.grid->mpi_communicator,&comparison);
  AssertThrowMPI(ierr);
  Assert (comparison == MPI_IDENT,ExcMessage("Matrix A and B must have a common MPI Communicator"));

  /*
   * The routine pgemr2d requires a BLACS context resembling at least the union of process grids
   * described by the BLACS contexts held by the ProcessGrids of matrix A and B.
   * As A and B share the same MPI communicator, there is no need to create a union MPI
   * communicator to initialise the BLACS context
   */
  int union_blacs_context = Csys2blacs_handle(this->grid->mpi_communicator);
  const char *order = "Col";
  int union_n_process_rows = Utilities::MPI::n_mpi_processes(this->grid->mpi_communicator);
  int union_n_process_columns = 1;
  Cblacs_gridinit(&union_blacs_context, order, union_n_process_rows, union_n_process_columns);

  int n_grid_rows_A,n_grid_columns_A,my_row_A,my_column_A;
  Cblacs_gridinfo(this->grid->blacs_context,&n_grid_rows_A,&n_grid_columns_A,&my_row_A,&my_column_A);

  //check whether process is in the BLACS context of matrix A
  const bool in_context_A = (my_row_A>=0 && my_row_A<n_grid_rows_A) &&
                            (my_column_A>=0 && my_column_A<n_grid_columns_A);


  int n_grid_rows_B,n_grid_columns_B,my_row_B,my_column_B;
  Cblacs_gridinfo(B.grid->blacs_context,&n_grid_rows_B,&n_grid_columns_B,&my_row_B,&my_column_B);

  //check whether process is in the BLACS context of matrix B
  const bool in_context_B = (my_row_B>=0 && my_row_B<n_grid_rows_B) &&
                            (my_column_B>=0 && my_column_B<n_grid_columns_B);

  const int n_rows_submatrix = submatrix_size.first;
  const int n_columns_submatrix =  submatrix_size.second;

  // due to Fortran indexing one has to be added
  int ia = offset_A.first+1, ja = offset_A.second+1;
  int ib = offset_B.first+1, jb = offset_B.second+1;

  std::array<int,9> desc_A, desc_B;

  const NumberType *loc_vals_A = nullptr;
  NumberType *loc_vals_B = nullptr;

  // Note: the function pgemr2d has to be called for all processes in the union BLACS context
  // If the calling process is not part of the BLACS context of A, desc_A[1] has to be -1
  // and all other parameters do not have to be set
  // If the calling process is not part of the BLACS context of B, desc_B[1] has to be -1
  // and all other parameters do not have to be set
  if (in_context_A)
    {
      if (this->values.size() != 0)
        loc_vals_A = & this->values[0];

      for (unsigned int i=0; i<desc_A.size(); ++i)
        desc_A[i] = this->descriptor[i];
    }
  else
    desc_A[1] =-1;

  if (in_context_B)
    {
      if (B.values.size() != 0)
        loc_vals_B = & B.values[0];

      for (unsigned int i=0; i<desc_B.size(); ++i)
        desc_B[i] = B.descriptor[i];
    }
  else
    desc_B[1]=-1;

  pgemr2d(&n_rows_submatrix, &n_columns_submatrix,
          loc_vals_A, &ia, &ja, desc_A.data(),
          loc_vals_B, &ib, &jb, desc_B.data(),
          &union_blacs_context);

  B.state = LAPACKSupport::matrix;

  //releasing the union BLACS context
  Cblacs_gridexit(union_blacs_context);
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
       * described by the BLACS contexts of matrix A and B
       */
      int union_blacs_context = Csys2blacs_handle(mpi_communicator_union);
      const char *order = "Col";
      int union_n_process_rows = Utilities::MPI::n_mpi_processes(mpi_communicator_union);
      int union_n_process_columns = 1;
      Cblacs_gridinit(&union_blacs_context, order, union_n_process_rows, union_n_process_columns);

      const NumberType *loc_vals_source = nullptr;
      NumberType *loc_vals_dest = nullptr;

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
void ScaLAPACKMatrix<NumberType>::add(const ScaLAPACKMatrix<NumberType> &B,
                                      const NumberType alpha,
                                      const NumberType beta,
                                      const bool transpose_B)
{
  if (transpose_B)
    {
      Assert (n_rows == B.n_columns, ExcDimensionMismatch(n_rows,B.n_columns));
      Assert (n_columns == B.n_rows, ExcDimensionMismatch(n_columns,B.n_rows));
      Assert(column_block_size==B.row_block_size,ExcDimensionMismatch(column_block_size,B.row_block_size));
      Assert(row_block_size==B.column_block_size,ExcDimensionMismatch(row_block_size,B.column_block_size));
    }
  else
    {
      Assert (n_rows == B.n_rows, ExcDimensionMismatch(n_rows,B.n_rows));
      Assert (n_columns == B.n_columns, ExcDimensionMismatch(n_columns,B.n_columns));
      Assert(column_block_size==B.column_block_size,ExcDimensionMismatch(column_block_size,B.column_block_size));
      Assert(row_block_size==B.row_block_size,ExcDimensionMismatch(row_block_size,B.row_block_size));
    }
  Assert(this->grid==B.grid,ExcMessage("The matrices A and B need to have the same process grid"));

  if (this->grid->mpi_process_is_active)
    {
      char trans_b = transpose_B ? 'T' : 'N';
      NumberType *A_loc = (this->values.size()>0) ? &this->values[0] : nullptr;
      const NumberType *B_loc = (B.values.size()>0) ? &B.values[0] : nullptr;

      pgeadd(&trans_b,&n_rows,&n_columns,
             &beta,B_loc,&B.submatrix_row,&B.submatrix_column,B.descriptor,
             &alpha,A_loc,&submatrix_row,&submatrix_column,descriptor);
    }
}



template <typename NumberType>
void add(const NumberType a,
		 const ScaLAPACKMatrix<NumberType> &B)
{
	add(B,1,a,false);
}



template <typename NumberType>
void Tadd(const NumberType a,
		  const ScaLAPACKMatrix<NumberType> &B)
{
	add(B,1,a,true);
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
                                                              std_cxx14::make_unique<ScaLAPACKMatrix<NumberType>>(n_rows,grid,row_block_size) :
                                                              std_cxx14::make_unique<ScaLAPACKMatrix<NumberType>>(grid->n_process_rows,grid->n_process_columns,grid,1,1);

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

  if (transpose)
    {
      Assert(n_columns==B.n_rows,ExcDimensionMismatch(n_columns,B.n_rows));
    }
  else
    {
      Assert(n_rows==B.n_rows,ExcDimensionMismatch(n_rows,B.n_rows));
    }

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

  if (property == LAPACKSupport::symmetric)
    return norm_symmetric(type);
  else
    return norm_general(type);
}



template <typename NumberType>
NumberType ScaLAPACKMatrix<NumberType>::linfty_norm() const
{
  const char type('I');

  if (property == LAPACKSupport::symmetric)
    return norm_symmetric(type);
  else
    return norm_general(type);
}



template <typename NumberType>
NumberType ScaLAPACKMatrix<NumberType>::frobenius_norm() const
{
  const char type('F');

  if (property == LAPACKSupport::symmetric)
    return norm_symmetric(type);
  else
    return norm_general(type);
}



template <typename NumberType>
NumberType ScaLAPACKMatrix<NumberType>::norm_general(const char type) const
{
  Assert (state == LAPACKSupport::matrix ||
          state == LAPACKSupport::inverse_matrix,
          ExcMessage("norms can be called in matrix state only."));
  Threads::Mutex::ScopedLock lock (mutex);
  NumberType res = 0.;

  if (grid->mpi_process_is_active)
    {
      const int iarow = indxg2p_(&submatrix_row, &row_block_size, &(grid->this_process_row), &first_process_row, &(grid->n_process_rows));
      const int iacol = indxg2p_(&submatrix_column, &column_block_size, &(grid->this_process_column), &first_process_column, &(grid->n_process_columns));
      const int mp0   = numroc_(&n_rows, &row_block_size, &(grid->this_process_row), &iarow, &(grid->n_process_rows));
      const int nq0   = numroc_(&n_columns, &column_block_size, &(grid->this_process_column), &iacol, &(grid->n_process_columns));

      // type='M': compute largest absolute value
      // type='F' || type='E': compute Frobenius norm
      // type='0' || type='1': compute infinity norm
      int lwork=0; // for type == 'M' || type == 'F' || type == 'E'
      if (type=='O' || type=='1')
        lwork = nq0;
      else if (type=='I')
        lwork = mp0;

      work.resize(lwork);
      const NumberType *A_loc = (this->values.size()>0) ? &this->values[0] : nullptr;
      res = plange(&type, &n_rows, &n_columns, A_loc, &submatrix_row, &submatrix_column, descriptor, work.data());
    }
  grid->send_to_inactive(&res);
  return res;
}



template <typename NumberType>
NumberType ScaLAPACKMatrix<NumberType>::norm_symmetric(const char type) const
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
      const NumberType *A_loc = (this->values.size()>0) ? &this->values[0] : nullptr;
      res = plansy(&type, &uplo, &n_columns, A_loc, &submatrix_row, &submatrix_column, descriptor, work.data());
    }
  grid->send_to_inactive(&res);
  return res;
}



template <typename NumberType>
void ScaLAPACKMatrix<NumberType>::save(const char *filename) const
{
#ifndef DEAL_II_WITH_HDF5
  (void)filename;
  AssertThrow(false, ExcMessage ("HDF5 support is disabled."));
#else
#  ifdef H5_HAVE_PARALLEL
  //implementation for configurations equipped with a parallel file system
  save_parallel(filename);

#  else
  //implementation for configurations with no parallel file system
  save_serial(filename);

#  endif
#endif
}



template <typename NumberType>
void ScaLAPACKMatrix<NumberType>::save_serial(const char *filename) const
{
#  ifndef DEAL_II_WITH_HDF5
  (void)filename;
  Assert(false,ExcInternalError());
#  else

  /*
   * The content of the distributed matrix is copied to a matrix using a 1x1 process grid.
   * Therefore, one process has all the data and can write it to a file.
   *
   * Create a 1x1 column grid which will be used to initialize
   * an effectively serial ScaLAPACK matrix to gather the contents from the current object
   */
  std::shared_ptr<Utilities::MPI::ProcessGrid> column_grid = std::make_shared<Utilities::MPI::ProcessGrid>(this->grid->mpi_communicator,1,1);

  const int MB=n_rows, NB=n_columns;
  ScaLAPACKMatrix<NumberType> tmp(n_rows,n_columns,column_grid,MB,NB);
  copy_to(tmp);

  // the 1x1 grid has only one process and this one writes
  // the content of the matrix to the HDF5 file
  if (tmp.grid->mpi_process_is_active)
    {
      herr_t status;

      // create a new file using default properties
      hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      // create the data space for the dataset
      hsize_t dims[2];
      //change order of rows and columns as ScaLAPACKMatrix uses column major ordering
      dims[0] = n_columns;
      dims[1] = n_rows;
      hid_t dataspace_id = H5Screate_simple(2, dims, nullptr);

      // create the dataset
      hid_t type_id = hdf5_type_id(&tmp.values[0]);
      hid_t dataset_id = H5Dcreate2(file_id, "/matrix",
                                    type_id, dataspace_id,
                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // write the dataset
      status = H5Dwrite(dataset_id, type_id,
                        H5S_ALL, H5S_ALL, H5P_DEFAULT,
                        &tmp.values[0]);
      AssertThrow(status >= 0, ExcIO());

      // end access to the dataset and release resources used by it
      status = H5Dclose(dataset_id);
      AssertThrow(status >= 0, ExcIO());

      // terminate access to the data space
      status = H5Sclose(dataspace_id);
      AssertThrow(status >= 0, ExcIO());

      // close the file.
      status = H5Fclose(file_id);
      AssertThrow(status >= 0, ExcIO());
    }
#  endif
}



template <typename NumberType>
void ScaLAPACKMatrix<NumberType>::save_parallel(const char *filename) const
{
#  ifndef DEAL_II_WITH_HDF5
  (void)filename;
  Assert(false,ExcInternalError());
#  else

  const unsigned int n_mpi_processes(Utilities::MPI::n_mpi_processes(this->grid->mpi_communicator));
  MPI_Info info = MPI_INFO_NULL;
  /*
   * The content of the distributed matrix is copied to a matrix using a 1xn_processes process grid.
   * Therefore, the processes hold contiguous chunks of the matrix, which they can write to the file
   *
   * Create a 1xn_processes column grid
  */
  std::shared_ptr<Utilities::MPI::ProcessGrid> column_grid = std::make_shared<Utilities::MPI::ProcessGrid>(this->grid->mpi_communicator,1,n_mpi_processes);

  const int MB=n_rows, NB=std::ceil(n_columns/n_mpi_processes);
  ScaLAPACKMatrix<NumberType> tmp(n_rows,n_columns,column_grid,MB,NB);
  copy_to(tmp);

  // get pointer to data held by the process
  NumberType *data = (tmp.values.size()>0) ? &tmp.values[0] : nullptr;

  herr_t status;
  // dataset dimensions
  hsize_t dims[2];

  // set up file access property list with parallel I/O access
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(plist_id, tmp.grid->mpi_communicator, info);
  AssertThrow(status >= 0, ExcIO());

  // create a new file collectively and release property list identifier
  hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  status = H5Pclose(plist_id);
  AssertThrow(status >= 0, ExcIO());

  // As ScaLAPACK, and therefore the class ScaLAPACKMatrix, uses column-major ordering
  // but HDF5 row-major ordering, we have to reverse entries related
  // to columns and rows in the following.
  // create the dataspace for the dataset
  dims[0] = tmp.n_columns;
  dims[1] = tmp.n_rows;

  hid_t filespace = H5Screate_simple(2, dims, nullptr);

  // create the dataset with default properties and close filespace
  hid_t type_id = hdf5_type_id(data);
  hid_t dset_id = H5Dcreate2(file_id, "/matrix", type_id,
                             filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Sclose(filespace);
  AssertThrow(status >= 0, ExcIO());

  // gather the number of local rows and columns from all processes
  std::vector<int> proc_n_local_rows(n_mpi_processes), proc_n_local_columns(n_mpi_processes);
  MPI_Allgather(&tmp.n_local_rows,1,MPI_INT,proc_n_local_rows.data(),1,MPI_INT,tmp.grid->mpi_communicator);
  MPI_Allgather(&tmp.n_local_columns,1,MPI_INT,proc_n_local_columns.data(),1,MPI_INT,tmp.grid->mpi_communicator);

  const unsigned int my_rank(Utilities::MPI::this_mpi_process(tmp.grid->mpi_communicator));

  // hyperslab selection parameters
  // each process defines dataset in memory and writes it to the hyperslab in the file
  hsize_t count[2];
  count[0] = tmp.n_local_columns;
  count[1] = tmp.n_rows;
  hid_t memspace = H5Screate_simple(2, count, nullptr);

  hsize_t offset[2] = {0};
  for (unsigned int i=0; i<my_rank; ++i)
    offset[0] += proc_n_local_columns[i];

  // select hyperslab in the file.
  filespace = H5Dget_space(dset_id);
  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, nullptr, count, nullptr);
  AssertThrow(status >= 0, ExcIO());

  // create property list for independent dataset write
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  AssertThrow(status >= 0, ExcIO());

  // process with no data will not participate in writing to the file
  if (tmp.values.size()>0)
    {
      status = H5Dwrite(dset_id, type_id, memspace, filespace,
                        plist_id, data);
      AssertThrow(status >= 0, ExcIO());
    }

  // close/release sources
  status = H5Dclose(dset_id);
  AssertThrow(status >= 0, ExcIO());
  status = H5Sclose(filespace);
  AssertThrow(status >= 0, ExcIO());
  status = H5Sclose(memspace);
  AssertThrow(status >= 0, ExcIO());
  status = H5Pclose(plist_id);
  AssertThrow(status >= 0, ExcIO());
  status = H5Fclose(file_id);
  AssertThrow(status >= 0, ExcIO());

#  endif
}



template <typename NumberType>
void ScaLAPACKMatrix<NumberType>::load(const char *filename)
{
#ifndef DEAL_II_WITH_HDF5
  (void)filename;
  AssertThrow(false, ExcMessage ("HDF5 support is disabled."));
#else
#  ifdef H5_HAVE_PARALLEL
  //implementation for configurations equipped with a parallel file system
  load_parallel(filename);
#  else
  //implementation for configurations with no parallel file system
  load_serial(filename);
#  endif
#endif
}



template <typename NumberType>
void ScaLAPACKMatrix<NumberType>::load_serial(const char *filename)
{
#  ifndef DEAL_II_WITH_HDF5
  (void)filename;
  Assert(false,ExcInternalError());
#  else

  /*
   * The content of the distributed matrix is copied to a matrix using a 1x1 process grid.
   * Therefore, one process has all the data and can write it to a file
   */
  //create a 1xP column grid with P being the number of MPI processes
  std::shared_ptr<Utilities::MPI::ProcessGrid> one_grid = std::make_shared<Utilities::MPI::ProcessGrid>(this->grid->mpi_communicator,1,1);

  const int MB=n_rows, NB=n_columns;
  ScaLAPACKMatrix<NumberType> tmp(n_rows,n_columns,one_grid,MB,NB);

  // the 1x1 grid has only one process and this one reads
  // the content of the matrix from the HDF5 file
  if (tmp.grid->mpi_process_is_active)
    {
      herr_t status;

      // open the file in read-only mode
      hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

      // open the dataset in the file
      hid_t dataset_id = H5Dopen2(file_id, "/matrix", H5P_DEFAULT);

      // check the datatype of the data in the file
      // datatype of source and destination must have the same class
      // see HDF User's Guide: 6.10. Data Transfer: Datatype Conversion and Selection
      hid_t datatype  = H5Dget_type(dataset_id);
      H5T_class_t t_class_in = H5Tget_class(datatype);
      H5T_class_t t_class = H5Tget_class(hdf5_type_id(&tmp.values[0]));
      AssertThrow(t_class_in == t_class,
                  ExcMessage("The data type of the matrix to be read does not match the archive"));

      // get dataspace handle
      hid_t dataspace_id = H5Dget_space(dataset_id);
      // get number of dimensions
      const int ndims = H5Sget_simple_extent_ndims(dataspace_id);
      AssertThrow(ndims==2, ExcIO());
      // get every dimension
      hsize_t dims[2];
      H5Sget_simple_extent_dims(dataspace_id, dims, nullptr);
      AssertThrow((int)dims[0]==n_columns,
                  ExcMessage("The number of columns of the matrix does not match the content of the archive"));
      AssertThrow((int)dims[1]==n_rows,
                  ExcMessage("The number of rows of the matrix does not match the content of the archive"));

      // read data
      status = H5Dread(dataset_id, hdf5_type_id(&tmp.values[0]), H5S_ALL, H5S_ALL,
                       H5P_DEFAULT, &tmp.values[0]);
      AssertThrow(status >= 0, ExcIO());

      // terminate access to the data space
      status = H5Sclose(dataspace_id);
      AssertThrow(status >= 0, ExcIO());

      // release data type handle
      status = H5Tclose(datatype);
      AssertThrow(status >= 0, ExcIO());

      // end access to the dataset and release resources used by it
      status = H5Dclose(dataset_id);
      AssertThrow(status >= 0, ExcIO());

      // close the file.
      status = H5Fclose(file_id);
      AssertThrow(status >= 0, ExcIO());
    }
  tmp.copy_to(*this);

#  endif // DEAL_II_WITH_HDF5
}



template <typename NumberType>
void ScaLAPACKMatrix<NumberType>::load_parallel(const char *filename)
{
#  ifndef DEAL_II_WITH_HDF5
  (void)filename;
  Assert(false,ExcInternalError());
#  else
#    ifndef H5_HAVE_PARALLEL
  Assert(false,ExcInternalError());
#    else

  const unsigned int n_mpi_processes(Utilities::MPI::n_mpi_processes(this->grid->mpi_communicator));
  MPI_Info info = MPI_INFO_NULL;
  /*
   * The content of the distributed matrix is copied to a matrix using a 1xn_processes process grid.
   * Therefore, the processes hold contiguous chunks of the matrix, which they can write to the file
   */
  //create a 1xP column grid with P being the number of MPI processes
  std::shared_ptr<Utilities::MPI::ProcessGrid> column_grid = std::make_shared<Utilities::MPI::ProcessGrid>(this->grid->mpi_communicator,1,n_mpi_processes);

  const int MB=n_rows, NB=std::ceil(n_columns/n_mpi_processes);
  ScaLAPACKMatrix<NumberType> tmp(n_rows,n_columns,column_grid,MB,NB);

  // get pointer to data held by the process
  NumberType *data = (tmp.values.size()>0) ? &tmp.values[0] : nullptr;

  herr_t status;

  // set up file access property list with parallel I/O access
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(plist_id, tmp.grid->mpi_communicator, info);
  AssertThrow(status >= 0, ExcIO());

  // open file collectively in read-only mode and release property list identifier
  hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, plist_id);
  status = H5Pclose(plist_id);
  AssertThrow(status >= 0, ExcIO());

  // open the dataset in the file collectively
  hid_t dataset_id = H5Dopen2(file_id, "/matrix", H5P_DEFAULT);

  // check the datatype of the dataset in the file
  // if the classes of type of the dataset and the matrix do not match abort
  // see HDF User's Guide: 6.10. Data Transfer: Datatype Conversion and Selection
  hid_t datatype = hdf5_type_id(data);
  hid_t datatype_inp  = H5Dget_type(dataset_id);
  H5T_class_t t_class_inp = H5Tget_class(datatype_inp);
  H5T_class_t t_class = H5Tget_class(datatype);
  AssertThrow(t_class_inp == t_class,
              ExcMessage("The data type of the matrix to be read does not match the archive"));

  // get the dimensions of the matrix stored in the file
  // get dataspace handle
  hid_t dataspace_id = H5Dget_space(dataset_id);
  // get number of dimensions
  const int ndims = H5Sget_simple_extent_ndims(dataspace_id);
  AssertThrow(ndims==2, ExcIO());
  // get every dimension
  hsize_t dims[2];
  status = H5Sget_simple_extent_dims(dataspace_id, dims, nullptr);
  AssertThrow(status >= 0, ExcIO());
  AssertThrow((int)dims[0]==n_columns,
              ExcMessage("The number of columns of the matrix does not match the content of the archive"));
  AssertThrow((int)dims[1]==n_rows,
              ExcMessage("The number of rows of the matrix does not match the content of the archive"));

  // gather the number of local rows and columns from all processes
  std::vector<int> proc_n_local_rows(n_mpi_processes), proc_n_local_columns(n_mpi_processes);
  MPI_Allgather(&tmp.n_local_rows,1,MPI_INT,proc_n_local_rows.data(),1,MPI_INT,tmp.grid->mpi_communicator);
  MPI_Allgather(&tmp.n_local_columns,1,MPI_INT,proc_n_local_columns.data(),1,MPI_INT,tmp.grid->mpi_communicator);

  const unsigned int my_rank(Utilities::MPI::this_mpi_process(tmp.grid->mpi_communicator));

  // hyperslab selection parameters
  // each process defines dataset in memory and writes it to the hyperslab in the file
  hsize_t count[2];
  count[0] = tmp.n_local_columns;
  count[1] = tmp.n_local_rows;

  hsize_t offset[2] = {0};
  for (unsigned int i=0; i<my_rank; ++i)
    offset[0] += proc_n_local_columns[i];

  // select hyperslab in the file
  status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
  AssertThrow(status >= 0, ExcIO());

  // create a memory dataspace independently
  hid_t memspace = H5Screate_simple(2, count, nullptr);

  // read data independently
  status = H5Dread(dataset_id, datatype, memspace, dataspace_id, H5P_DEFAULT, data);
  AssertThrow(status >= 0, ExcIO());

  // close/release sources
  status = H5Sclose(memspace);
  AssertThrow(status >= 0, ExcIO());
  status = H5Dclose(dataset_id);
  AssertThrow(status >= 0, ExcIO());
  status = H5Sclose(dataspace_id);
  AssertThrow(status >= 0, ExcIO());
  status = H5Fclose(file_id);
  AssertThrow(status >= 0, ExcIO());

  // copying the distributed matrices
  tmp.copy_to(*this);

#    endif // H5_HAVE_PARALLEL
#  endif // DEAL_II_WITH_HDF5
}



template <typename NumberType>
void ScaLAPACKMatrix<NumberType>::scale_columns(std::vector<NumberType> &factors)
{
  Assert(n_columns==(int)factors.size(),ExcDimensionMismatch(n_columns,factors.size()));

  if (this->grid->mpi_process_is_active)
    {
      for (int i=0; i<n_local_rows; ++i)
        {
          for (int j=0; j<n_local_columns; ++j)
            {
              const int glob_j = global_column(j);
              local_el(i,j) *= factors[glob_j];
            }
        }
    }
}



template <typename NumberType>
void ScaLAPACKMatrix<NumberType>::scale_rows(std::vector<NumberType> &factors)
{
  Assert(n_rows==(int)factors.size(),ExcDimensionMismatch(n_rows,factors.size()));

  if (this->grid->mpi_process_is_active)
    {
      for (int i=0; i<n_local_rows; ++i)
        {
          const int glob_i = global_row(i);
          for (int j=0; j<n_local_columns; ++j)
            {
              local_el(i,j) *= factors[glob_i];
            }
        }
    }
}



// instantiations
template class ScaLAPACKMatrix<double>;
template class ScaLAPACKMatrix<float>;


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SCALAPACK
