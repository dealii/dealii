// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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


#include <deal.II/lac/scalapack.h>

#ifdef DEAL_II_WITH_SCALAPACK

#  include <deal.II/base/array_view.h>
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/mpi.templates.h>

#  include <deal.II/lac/scalapack.templates.h>

#  ifdef DEAL_II_WITH_HDF5
#    include <hdf5.h>
#  endif

#  include <memory>

DEAL_II_NAMESPACE_OPEN

#  ifdef DEAL_II_WITH_HDF5

template <typename number>
inline hid_t
hdf5_type_id(const number *)
{
  Assert(false, dealii::ExcNotImplemented());
  // don't know what to put here; it does not matter
  return -1;
}

inline hid_t
hdf5_type_id(const double *)
{
  return H5T_NATIVE_DOUBLE;
}

inline hid_t
hdf5_type_id(const float *)
{
  return H5T_NATIVE_FLOAT;
}

inline hid_t
hdf5_type_id(const int *)
{
  return H5T_NATIVE_INT;
}

inline hid_t
hdf5_type_id(const unsigned int *)
{
  return H5T_NATIVE_UINT;
}

inline hid_t
hdf5_type_id(const char *)
{
  return H5T_NATIVE_CHAR;
}
#  endif // DEAL_II_WITH_HDF5



template <typename NumberType>
ScaLAPACKMatrix<NumberType>::ScaLAPACKMatrix(
  const size_type                                           n_rows_,
  const size_type                                           n_columns_,
  const std::shared_ptr<const Utilities::MPI::ProcessGrid> &process_grid,
  const size_type                                           row_block_size_,
  const size_type                                           column_block_size_,
  const LAPACKSupport::Property                             property_)
  : uplo('L')
  , // for non-symmetric matrices this is not needed
  first_process_row(0)
  , first_process_column(0)
  , submatrix_row(1)
  , submatrix_column(1)
{
  reinit(n_rows_,
         n_columns_,
         process_grid,
         row_block_size_,
         column_block_size_,
         property_);
}



template <typename NumberType>
ScaLAPACKMatrix<NumberType>::ScaLAPACKMatrix(
  const size_type                                           size,
  const std::shared_ptr<const Utilities::MPI::ProcessGrid> &process_grid,
  const size_type                                           block_size,
  const LAPACKSupport::Property                             property)
  : ScaLAPACKMatrix<NumberType>(size,
                                size,
                                process_grid,
                                block_size,
                                block_size,
                                property)
{}



template <typename NumberType>
ScaLAPACKMatrix<NumberType>::ScaLAPACKMatrix(
  const std::string &                                       filename,
  const std::shared_ptr<const Utilities::MPI::ProcessGrid> &process_grid,
  const size_type                                           row_block_size,
  const size_type                                           column_block_size)
  : uplo('L')
  , // for non-symmetric matrices this is not needed
  first_process_row(0)
  , first_process_column(0)
  , submatrix_row(1)
  , submatrix_column(1)
{
#  ifndef DEAL_II_WITH_HDF5
  (void)filename;
  (void)process_grid;
  (void)row_block_size;
  (void)column_block_size;
  Assert(
    false,
    ExcMessage(
      "This function is only available when deal.II is configured with HDF5"));
#  else

  const unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(process_grid->mpi_communicator));

  // Before reading the content from disk the root process determines the
  // dimensions of the matrix. Subsequently, memory is allocated by a call to
  // reinit() and the matrix is loaded by a call to load().
  if (this_mpi_process == 0)
    {
      herr_t status = 0;

      // open file in read-only mode
      hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      AssertThrow(file >= 0, ExcIO());

      // get data set in file
      hid_t dataset = H5Dopen2(file, "/matrix", H5P_DEFAULT);
      AssertThrow(dataset >= 0, ExcIO());

      // determine file space
      hid_t filespace = H5Dget_space(dataset);

      // get number of dimensions in data set
      int rank = H5Sget_simple_extent_ndims(filespace);
      AssertThrow(rank == 2, ExcIO());
      hsize_t dims[2];
      status = H5Sget_simple_extent_dims(filespace, dims, nullptr);
      AssertThrow(status >= 0, ExcIO());

      // due to ScaLAPACK's column-major memory layout the dimensions are
      // swapped
      n_rows    = dims[1];
      n_columns = dims[0];

      // close/release resources
      status = H5Sclose(filespace);
      AssertThrow(status >= 0, ExcIO());
      status = H5Dclose(dataset);
      AssertThrow(status >= 0, ExcIO());
      status = H5Fclose(file);
      AssertThrow(status >= 0, ExcIO());
    }
  int ierr = MPI_Bcast(&n_rows,
                       1,
                       Utilities::MPI::internal::mpi_type_id(&n_rows),
                       0 /*from root*/,
                       process_grid->mpi_communicator);
  AssertThrowMPI(ierr);

  ierr = MPI_Bcast(&n_columns,
                   1,
                   Utilities::MPI::internal::mpi_type_id(&n_columns),
                   0 /*from root*/,
                   process_grid->mpi_communicator);
  AssertThrowMPI(ierr);

  // the property will be overwritten by the subsequent call to load()
  reinit(n_rows,
         n_columns,
         process_grid,
         row_block_size,
         column_block_size,
         LAPACKSupport::Property::general);

  load(filename.c_str());

#  endif // DEAL_II_WITH_HDF5
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::reinit(
  const size_type                                           n_rows_,
  const size_type                                           n_columns_,
  const std::shared_ptr<const Utilities::MPI::ProcessGrid> &process_grid,
  const size_type                                           row_block_size_,
  const size_type                                           column_block_size_,
  const LAPACKSupport::Property                             property_)
{
  Assert(row_block_size_ > 0, ExcMessage("Row block size has to be positive."));
  Assert(column_block_size_ > 0,
         ExcMessage("Column block size has to be positive."));
  Assert(
    row_block_size_ <= n_rows_,
    ExcMessage(
      "Row block size can not be greater than the number of rows of the matrix"));
  Assert(
    column_block_size_ <= n_columns_,
    ExcMessage(
      "Column block size can not be greater than the number of columns of the matrix"));

  state             = LAPACKSupport::State::matrix;
  property          = property_;
  grid              = process_grid;
  n_rows            = n_rows_;
  n_columns         = n_columns_;
  row_block_size    = row_block_size_;
  column_block_size = column_block_size_;

  if (grid->mpi_process_is_active)
    {
      // Get local sizes:
      n_local_rows    = numroc_(&n_rows,
                             &row_block_size,
                             &(grid->this_process_row),
                             &first_process_row,
                             &(grid->n_process_rows));
      n_local_columns = numroc_(&n_columns,
                                &column_block_size,
                                &(grid->this_process_column),
                                &first_process_column,
                                &(grid->n_process_columns));

      // LLD_A = MAX(1,NUMROC(M_A, MB_A, MYROW, RSRC_A, NPROW)), different
      // between processes
      int lda = std::max(1, n_local_rows);

      int info = 0;
      descinit_(descriptor,
                &n_rows,
                &n_columns,
                &row_block_size,
                &column_block_size,
                &first_process_row,
                &first_process_column,
                &(grid->blacs_context),
                &lda,
                &info);
      AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("descinit_", info));

      this->TransposeTable<NumberType>::reinit(n_local_rows, n_local_columns);
    }
  else
    {
      // set process-local variables to something telling:
      n_local_rows    = -1;
      n_local_columns = -1;
      std::fill(std::begin(descriptor), std::end(descriptor), -1);
    }
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::reinit(
  const size_type                                           size,
  const std::shared_ptr<const Utilities::MPI::ProcessGrid> &process_grid,
  const size_type                                           block_size,
  const LAPACKSupport::Property                             property)
{
  reinit(size, size, process_grid, block_size, block_size, property);
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::set_property(
  const LAPACKSupport::Property property_)
{
  property = property_;
}



template <typename NumberType>
LAPACKSupport::Property
ScaLAPACKMatrix<NumberType>::get_property() const
{
  return property;
}



template <typename NumberType>
LAPACKSupport::State
ScaLAPACKMatrix<NumberType>::get_state() const
{
  return state;
}



template <typename NumberType>
ScaLAPACKMatrix<NumberType> &
ScaLAPACKMatrix<NumberType>::operator=(const FullMatrix<NumberType> &matrix)
{
  // FIXME: another way to copy is to use pdgeadd_ PBLAS routine.
  // This routine computes the sum of two matrices B := a*A + b*B.
  // Matrices can have different distribution,in particular matrix A can
  // be owned by only one process, so we can set a=1 and b=0 to copy
  // non-distributed matrix A into distributed matrix B.
  Assert(n_rows == int(matrix.m()), ExcDimensionMismatch(n_rows, matrix.m()));
  Assert(n_columns == int(matrix.n()),
         ExcDimensionMismatch(n_columns, matrix.n()));

  if (grid->mpi_process_is_active)
    {
      for (int i = 0; i < n_local_rows; ++i)
        {
          const int glob_i = global_row(i);
          for (int j = 0; j < n_local_columns; ++j)
            {
              const int glob_j = global_column(j);
              local_el(i, j)   = matrix(glob_i, glob_j);
            }
        }
    }
  state = LAPACKSupport::matrix;
  return *this;
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::copy_from(const LAPACKFullMatrix<NumberType> &B,
                                       const unsigned int                  rank)
{
  if (n_rows * n_columns == 0)
    return;

  const unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(this->grid->mpi_communicator));

#  ifdef DEBUG
  Assert(Utilities::MPI::max(rank, this->grid->mpi_communicator) == rank,
         ExcMessage("All processes have to call routine with identical rank"));
  Assert(Utilities::MPI::min(rank, this->grid->mpi_communicator) == rank,
         ExcMessage("All processes have to call routine with identical rank"));
#  endif

  // root process has to be active in the grid of A
  if (this_mpi_process == rank)
    {
      Assert(grid->mpi_process_is_active, ExcInternalError());
      Assert(n_rows == int(B.m()), ExcDimensionMismatch(n_rows, B.m()));
      Assert(n_columns == int(B.n()), ExcDimensionMismatch(n_columns, B.n()));
    }
  // Create 1x1 grid for matrix B.
  // The underlying grid for matrix B only contains the process #rank.
  // This grid will be used to copy the serial matrix B to the distributed
  // matrix using the ScaLAPACK routine pgemr2d.
  MPI_Group group_A;
  MPI_Comm_group(this->grid->mpi_communicator, &group_A);
  const int              n = 1;
  const std::vector<int> ranks(n, rank);
  MPI_Group              group_B;
  MPI_Group_incl(group_A, n, DEAL_II_MPI_CONST_CAST(ranks.data()), &group_B);
  MPI_Comm communicator_B;

  const int mpi_tag = Utilities::MPI::internal::Tags::scalapack_copy_from;
  Utilities::MPI::create_group(this->grid->mpi_communicator,
                               group_B,
                               mpi_tag,
                               &communicator_B);
  int n_proc_rows_B = 1, n_proc_cols_B = 1;
  int this_process_row_B = -1, this_process_column_B = -1;
  int blacs_context_B = -1;
  if (MPI_COMM_NULL != communicator_B)
    {
      // Initialize Cblas context from the provided communicator
      blacs_context_B   = Csys2blacs_handle(communicator_B);
      const char *order = "Col";
      Cblacs_gridinit(&blacs_context_B, order, n_proc_rows_B, n_proc_cols_B);
      Cblacs_gridinfo(blacs_context_B,
                      &n_proc_rows_B,
                      &n_proc_cols_B,
                      &this_process_row_B,
                      &this_process_column_B);
      Assert(n_proc_rows_B * n_proc_cols_B == 1, ExcInternalError());
      // the active process of grid B has to be process #rank of the
      // communicator attached to A
      Assert(this_mpi_process == rank, ExcInternalError());
    }
  const bool mpi_process_is_active_B =
    (this_process_row_B >= 0 && this_process_column_B >= 0);

  // create descriptor for matrix B
  std::vector<int> descriptor_B(9, -1);
  const int        first_process_row_B = 0, first_process_col_B = 0;

  if (mpi_process_is_active_B)
    {
      // Get local sizes:
      int n_local_rows_B = numroc_(&n_rows,
                                   &n_rows,
                                   &this_process_row_B,
                                   &first_process_row_B,
                                   &n_proc_rows_B);
      int n_local_cols_B = numroc_(&n_columns,
                                   &n_columns,
                                   &this_process_column_B,
                                   &first_process_col_B,
                                   &n_proc_cols_B);
      Assert(n_local_rows_B == n_rows, ExcInternalError());
      Assert(n_local_cols_B == n_columns, ExcInternalError());
      (void)n_local_cols_B;

      int lda  = std::max(1, n_local_rows_B);
      int info = 0;
      descinit_(descriptor_B.data(),
                &n_rows,
                &n_columns,
                &n_rows,
                &n_columns,
                &first_process_row_B,
                &first_process_col_B,
                &blacs_context_B,
                &lda,
                &info);
      AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("descinit_", info));
    }
  if (this->grid->mpi_process_is_active)
    {
      const int   ii = 1;
      NumberType *loc_vals_A =
        this->values.size() > 0 ? this->values.data() : nullptr;
      const NumberType *loc_vals_B =
        mpi_process_is_active_B ? &(B(0, 0)) : nullptr;

      // pgemr2d has to be called only for processes active on grid attached to
      // matrix A
      pgemr2d(&n_rows,
              &n_columns,
              loc_vals_B,
              &ii,
              &ii,
              descriptor_B.data(),
              loc_vals_A,
              &ii,
              &ii,
              this->descriptor,
              &(this->grid->blacs_context));
    }
  if (mpi_process_is_active_B)
    Cblacs_gridexit(blacs_context_B);

  MPI_Group_free(&group_A);
  MPI_Group_free(&group_B);
  if (MPI_COMM_NULL != communicator_B)
    MPI_Comm_free(&communicator_B);

  state = LAPACKSupport::matrix;
}



template <typename NumberType>
unsigned int
ScaLAPACKMatrix<NumberType>::global_row(const unsigned int loc_row) const
{
  Assert(n_local_rows >= 0 && loc_row < static_cast<unsigned int>(n_local_rows),
         ExcIndexRange(loc_row, 0, n_local_rows));
  const int i = loc_row + 1;
  return indxl2g_(&i,
                  &row_block_size,
                  &(grid->this_process_row),
                  &first_process_row,
                  &(grid->n_process_rows)) -
         1;
}



template <typename NumberType>
unsigned int
ScaLAPACKMatrix<NumberType>::global_column(const unsigned int loc_column) const
{
  Assert(n_local_columns >= 0 &&
           loc_column < static_cast<unsigned int>(n_local_columns),
         ExcIndexRange(loc_column, 0, n_local_columns));
  const int j = loc_column + 1;
  return indxl2g_(&j,
                  &column_block_size,
                  &(grid->this_process_column),
                  &first_process_column,
                  &(grid->n_process_columns)) -
         1;
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::copy_to(LAPACKFullMatrix<NumberType> &B,
                                     const unsigned int            rank) const
{
  if (n_rows * n_columns == 0)
    return;

  const unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(this->grid->mpi_communicator));

#  ifdef DEBUG
  Assert(Utilities::MPI::max(rank, this->grid->mpi_communicator) == rank,
         ExcMessage("All processes have to call routine with identical rank"));
  Assert(Utilities::MPI::min(rank, this->grid->mpi_communicator) == rank,
         ExcMessage("All processes have to call routine with identical rank"));
#  endif

  if (this_mpi_process == rank)
    {
      // the process which gets the serial copy has to be in the process grid
      Assert(this->grid->is_process_active(), ExcInternalError());
      Assert(n_rows == int(B.m()), ExcDimensionMismatch(n_rows, B.m()));
      Assert(n_columns == int(B.n()), ExcDimensionMismatch(n_columns, B.n()));
    }

  // Create 1x1 grid for matrix B.
  // The underlying grid for matrix B only contains the process #rank.
  // This grid will be used to copy to the distributed matrix to the serial
  // matrix B using the ScaLAPACK routine pgemr2d.
  MPI_Group group_A;
  MPI_Comm_group(this->grid->mpi_communicator, &group_A);
  const int              n = 1;
  const std::vector<int> ranks(n, rank);
  MPI_Group              group_B;
  MPI_Group_incl(group_A, n, DEAL_II_MPI_CONST_CAST(ranks.data()), &group_B);
  MPI_Comm communicator_B;

  const int mpi_tag = Utilities::MPI::internal::Tags::scalapack_copy_to;
  Utilities::MPI::create_group(this->grid->mpi_communicator,
                               group_B,
                               mpi_tag,
                               &communicator_B);
  int n_proc_rows_B = 1, n_proc_cols_B = 1;
  int this_process_row_B = -1, this_process_column_B = -1;
  int blacs_context_B = -1;
  if (MPI_COMM_NULL != communicator_B)
    {
      // Initialize Cblas context from the provided communicator
      blacs_context_B   = Csys2blacs_handle(communicator_B);
      const char *order = "Col";
      Cblacs_gridinit(&blacs_context_B, order, n_proc_rows_B, n_proc_cols_B);
      Cblacs_gridinfo(blacs_context_B,
                      &n_proc_rows_B,
                      &n_proc_cols_B,
                      &this_process_row_B,
                      &this_process_column_B);
      Assert(n_proc_rows_B * n_proc_cols_B == 1, ExcInternalError());
      // the active process of grid B has to be process #rank of the
      // communicator attached to A
      Assert(this_mpi_process == rank, ExcInternalError());
    }
  const bool mpi_process_is_active_B =
    (this_process_row_B >= 0 && this_process_column_B >= 0);

  // create descriptor for matrix B
  std::vector<int> descriptor_B(9, -1);
  const int        first_process_row_B = 0, first_process_col_B = 0;

  if (mpi_process_is_active_B)
    {
      // Get local sizes:
      int n_local_rows_B = numroc_(&n_rows,
                                   &n_rows,
                                   &this_process_row_B,
                                   &first_process_row_B,
                                   &n_proc_rows_B);
      int n_local_cols_B = numroc_(&n_columns,
                                   &n_columns,
                                   &this_process_column_B,
                                   &first_process_col_B,
                                   &n_proc_cols_B);
      Assert(n_local_rows_B == n_rows, ExcInternalError());
      Assert(n_local_cols_B == n_columns, ExcInternalError());
      (void)n_local_cols_B;

      int lda  = std::max(1, n_local_rows_B);
      int info = 0;
      // fill descriptor for matrix B
      descinit_(descriptor_B.data(),
                &n_rows,
                &n_columns,
                &n_rows,
                &n_columns,
                &first_process_row_B,
                &first_process_col_B,
                &blacs_context_B,
                &lda,
                &info);
      AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("descinit_", info));
    }
  // pgemr2d has to be called only for processes active on grid attached to
  // matrix A
  if (this->grid->mpi_process_is_active)
    {
      const int         ii = 1;
      const NumberType *loc_vals_A =
        this->values.size() > 0 ? this->values.data() : nullptr;
      NumberType *loc_vals_B = mpi_process_is_active_B ? &(B(0, 0)) : nullptr;

      pgemr2d(&n_rows,
              &n_columns,
              loc_vals_A,
              &ii,
              &ii,
              this->descriptor,
              loc_vals_B,
              &ii,
              &ii,
              descriptor_B.data(),
              &(this->grid->blacs_context));
    }
  if (mpi_process_is_active_B)
    Cblacs_gridexit(blacs_context_B);

  MPI_Group_free(&group_A);
  MPI_Group_free(&group_B);
  if (MPI_COMM_NULL != communicator_B)
    MPI_Comm_free(&communicator_B);
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::copy_to(FullMatrix<NumberType> &matrix) const
{
  // FIXME: use PDGEMR2D for copying?
  // PDGEMR2D copies a submatrix of A on a submatrix of B.
  // A and B can have different distributions
  // see http://icl.cs.utk.edu/lapack-forum/viewtopic.php?t=50
  Assert(n_rows == int(matrix.m()), ExcDimensionMismatch(n_rows, matrix.m()));
  Assert(n_columns == int(matrix.n()),
         ExcDimensionMismatch(n_columns, matrix.n()));

  matrix = 0.;
  if (grid->mpi_process_is_active)
    {
      for (int i = 0; i < n_local_rows; ++i)
        {
          const int glob_i = global_row(i);
          for (int j = 0; j < n_local_columns; ++j)
            {
              const int glob_j       = global_column(j);
              matrix(glob_i, glob_j) = local_el(i, j);
            }
        }
    }
  Utilities::MPI::sum(matrix, grid->mpi_communicator, matrix);

  // we could move the following lines under the main loop above,
  // but they would be dependent on glob_i and glob_j, which
  // won't make it much prettier
  if (state == LAPACKSupport::cholesky)
    {
      if (property == LAPACKSupport::lower_triangular)
        for (unsigned int i = 0; i < matrix.n(); ++i)
          for (unsigned int j = i + 1; j < matrix.m(); ++j)
            matrix(i, j) = 0.;
      else if (property == LAPACKSupport::upper_triangular)
        for (unsigned int i = 0; i < matrix.n(); ++i)
          for (unsigned int j = 0; j < i; ++j)
            matrix(i, j) = 0.;
    }
  else if (property == LAPACKSupport::symmetric &&
           state == LAPACKSupport::inverse_matrix)
    {
      if (uplo == 'L')
        for (unsigned int i = 0; i < matrix.n(); ++i)
          for (unsigned int j = i + 1; j < matrix.m(); ++j)
            matrix(i, j) = matrix(j, i);
      else if (uplo == 'U')
        for (unsigned int i = 0; i < matrix.n(); ++i)
          for (unsigned int j = 0; j < i; ++j)
            matrix(i, j) = matrix(j, i);
    }
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::copy_to(
  ScaLAPACKMatrix<NumberType> &                B,
  const std::pair<unsigned int, unsigned int> &offset_A,
  const std::pair<unsigned int, unsigned int> &offset_B,
  const std::pair<unsigned int, unsigned int> &submatrix_size) const
{
  // submatrix is empty
  if (submatrix_size.first == 0 || submatrix_size.second == 0)
    return;

  // range checking for matrix A
  AssertIndexRange(offset_A.first, n_rows - submatrix_size.first + 1);
  AssertIndexRange(offset_A.second, n_columns - submatrix_size.second + 1);

  // range checking for matrix B
  AssertIndexRange(offset_B.first, B.n_rows - submatrix_size.first + 1);
  AssertIndexRange(offset_B.second, B.n_columns - submatrix_size.second + 1);

  // Currently, copying of matrices will only be supported if A and B share the
  // same MPI communicator
  int ierr, comparison;
  ierr = MPI_Comm_compare(grid->mpi_communicator,
                          B.grid->mpi_communicator,
                          &comparison);
  AssertThrowMPI(ierr);
  Assert(comparison == MPI_IDENT,
         ExcMessage("Matrix A and B must have a common MPI Communicator"));

  /*
   * The routine pgemr2d requires a BLACS context resembling at least the union
   * of process grids described by the BLACS contexts held by the ProcessGrids
   * of matrix A and B. As A and B share the same MPI communicator, there is no
   * need to create a union MPI communicator to initialise the BLACS context
   */
  int union_blacs_context = Csys2blacs_handle(this->grid->mpi_communicator);
  const char *order       = "Col";
  int         union_n_process_rows =
    Utilities::MPI::n_mpi_processes(this->grid->mpi_communicator);
  int union_n_process_columns = 1;
  Cblacs_gridinit(&union_blacs_context,
                  order,
                  union_n_process_rows,
                  union_n_process_columns);

  int n_grid_rows_A, n_grid_columns_A, my_row_A, my_column_A;
  Cblacs_gridinfo(this->grid->blacs_context,
                  &n_grid_rows_A,
                  &n_grid_columns_A,
                  &my_row_A,
                  &my_column_A);

  // check whether process is in the BLACS context of matrix A
  const bool in_context_A =
    (my_row_A >= 0 && my_row_A < n_grid_rows_A) &&
    (my_column_A >= 0 && my_column_A < n_grid_columns_A);

  int n_grid_rows_B, n_grid_columns_B, my_row_B, my_column_B;
  Cblacs_gridinfo(B.grid->blacs_context,
                  &n_grid_rows_B,
                  &n_grid_columns_B,
                  &my_row_B,
                  &my_column_B);

  // check whether process is in the BLACS context of matrix B
  const bool in_context_B =
    (my_row_B >= 0 && my_row_B < n_grid_rows_B) &&
    (my_column_B >= 0 && my_column_B < n_grid_columns_B);

  const int n_rows_submatrix    = submatrix_size.first;
  const int n_columns_submatrix = submatrix_size.second;

  // due to Fortran indexing one has to be added
  int ia = offset_A.first + 1, ja = offset_A.second + 1;
  int ib = offset_B.first + 1, jb = offset_B.second + 1;

  std::array<int, 9> desc_A, desc_B;

  const NumberType *loc_vals_A = nullptr;
  NumberType *      loc_vals_B = nullptr;

  // Note: the function pgemr2d has to be called for all processes in the union
  // BLACS context If the calling process is not part of the BLACS context of A,
  // desc_A[1] has to be -1 and all other parameters do not have to be set If
  // the calling process is not part of the BLACS context of B, desc_B[1] has to
  // be -1 and all other parameters do not have to be set
  if (in_context_A)
    {
      if (this->values.size() != 0)
        loc_vals_A = this->values.data();

      for (unsigned int i = 0; i < desc_A.size(); ++i)
        desc_A[i] = this->descriptor[i];
    }
  else
    desc_A[1] = -1;

  if (in_context_B)
    {
      if (B.values.size() != 0)
        loc_vals_B = B.values.data();

      for (unsigned int i = 0; i < desc_B.size(); ++i)
        desc_B[i] = B.descriptor[i];
    }
  else
    desc_B[1] = -1;

  pgemr2d(&n_rows_submatrix,
          &n_columns_submatrix,
          loc_vals_A,
          &ia,
          &ja,
          desc_A.data(),
          loc_vals_B,
          &ib,
          &jb,
          desc_B.data(),
          &union_blacs_context);

  B.state = LAPACKSupport::matrix;

  // releasing the union BLACS context
  Cblacs_gridexit(union_blacs_context);
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::copy_to(ScaLAPACKMatrix<NumberType> &dest) const
{
  Assert(n_rows == dest.n_rows, ExcDimensionMismatch(n_rows, dest.n_rows));
  Assert(n_columns == dest.n_columns,
         ExcDimensionMismatch(n_columns, dest.n_columns));

  if (this->grid->mpi_process_is_active)
    AssertThrow(
      this->descriptor[0] == 1,
      ExcMessage(
        "Copying of ScaLAPACK matrices only implemented for dense matrices"));
  if (dest.grid->mpi_process_is_active)
    AssertThrow(
      dest.descriptor[0] == 1,
      ExcMessage(
        "Copying of ScaLAPACK matrices only implemented for dense matrices"));

  /*
   * just in case of different process grids or block-cyclic distributions
   * inter-process communication is necessary
   * if distributed matrices have the same process grid and block sizes, local
   * copying is enough
   */
  if ((this->grid != dest.grid) || (row_block_size != dest.row_block_size) ||
      (column_block_size != dest.column_block_size))
    {
      /*
       * get the MPI communicator, which is the union of the source and
       * destination MPI communicator
       */
      int       ierr = 0;
      MPI_Group group_source, group_dest, group_union;
      ierr = MPI_Comm_group(this->grid->mpi_communicator, &group_source);
      AssertThrowMPI(ierr);
      ierr = MPI_Comm_group(dest.grid->mpi_communicator, &group_dest);
      AssertThrowMPI(ierr);
      ierr = MPI_Group_union(group_source, group_dest, &group_union);
      AssertThrowMPI(ierr);
      MPI_Comm mpi_communicator_union;

      // to create a communicator representing the union of the source
      // and destination MPI communicator we need a communicator that
      // is guaranteed to contain all desired processes -- i.e.,
      // MPI_COMM_WORLD. on the other hand, as documented in the MPI
      // standard, MPI_Comm_create_group is not collective on all
      // processes in the first argument, but instead is collective on
      // only those processes listed in the group. in other words,
      // there is really no harm in passing MPI_COMM_WORLD as the
      // first argument, even if the program we are currently running
      // and that is calling this function only works on a subset of
      // processes. the same holds for the wrapper/fallback we are using here.

      const int mpi_tag = Utilities::MPI::internal::Tags::scalapack_copy_to2;
      ierr              = Utilities::MPI::create_group(MPI_COMM_WORLD,
                                          group_union,
                                          mpi_tag,
                                          &mpi_communicator_union);
      AssertThrowMPI(ierr);

      /*
       * The routine pgemr2d requires a BLACS context resembling at least the
       * union of process grids described by the BLACS contexts of matrix A and
       * B
       */
      int union_blacs_context = Csys2blacs_handle(mpi_communicator_union);
      const char *order       = "Col";
      int         union_n_process_rows =
        Utilities::MPI::n_mpi_processes(mpi_communicator_union);
      int union_n_process_columns = 1;
      Cblacs_gridinit(&union_blacs_context,
                      order,
                      union_n_process_rows,
                      union_n_process_columns);

      const NumberType *loc_vals_source = nullptr;
      NumberType *      loc_vals_dest   = nullptr;

      if (this->grid->mpi_process_is_active && (this->values.size() > 0))
        {
          AssertThrow(this->values.size() > 0,
                      dealii::ExcMessage(
                        "source: process is active but local matrix empty"));
          loc_vals_source = this->values.data();
        }
      if (dest.grid->mpi_process_is_active && (dest.values.size() > 0))
        {
          AssertThrow(
            dest.values.size() > 0,
            dealii::ExcMessage(
              "destination: process is active but local matrix empty"));
          loc_vals_dest = dest.values.data();
        }
      pgemr2d(&n_rows,
              &n_columns,
              loc_vals_source,
              &submatrix_row,
              &submatrix_column,
              descriptor,
              loc_vals_dest,
              &dest.submatrix_row,
              &dest.submatrix_column,
              dest.descriptor,
              &union_blacs_context);

      Cblacs_gridexit(union_blacs_context);

      if (mpi_communicator_union != MPI_COMM_NULL)
        {
          ierr = MPI_Comm_free(&mpi_communicator_union);
          AssertThrowMPI(ierr);
        }
      ierr = MPI_Group_free(&group_source);
      AssertThrowMPI(ierr);
      ierr = MPI_Group_free(&group_dest);
      AssertThrowMPI(ierr);
      ierr = MPI_Group_free(&group_union);
      AssertThrowMPI(ierr);
    }
  else
    // process is active in the process grid
    if (this->grid->mpi_process_is_active)
    dest.values = this->values;

  dest.state    = state;
  dest.property = property;
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::copy_transposed(
  const ScaLAPACKMatrix<NumberType> &B)
{
  add(B, 0, 1, true);
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::add(const ScaLAPACKMatrix<NumberType> &B,
                                 const NumberType                   alpha,
                                 const NumberType                   beta,
                                 const bool                         transpose_B)
{
  if (transpose_B)
    {
      Assert(n_rows == B.n_columns, ExcDimensionMismatch(n_rows, B.n_columns));
      Assert(n_columns == B.n_rows, ExcDimensionMismatch(n_columns, B.n_rows));
      Assert(column_block_size == B.row_block_size,
             ExcDimensionMismatch(column_block_size, B.row_block_size));
      Assert(row_block_size == B.column_block_size,
             ExcDimensionMismatch(row_block_size, B.column_block_size));
    }
  else
    {
      Assert(n_rows == B.n_rows, ExcDimensionMismatch(n_rows, B.n_rows));
      Assert(n_columns == B.n_columns,
             ExcDimensionMismatch(n_columns, B.n_columns));
      Assert(column_block_size == B.column_block_size,
             ExcDimensionMismatch(column_block_size, B.column_block_size));
      Assert(row_block_size == B.row_block_size,
             ExcDimensionMismatch(row_block_size, B.row_block_size));
    }
  Assert(this->grid == B.grid,
         ExcMessage("The matrices A and B need to have the same process grid"));

  if (this->grid->mpi_process_is_active)
    {
      char        trans_b = transpose_B ? 'T' : 'N';
      NumberType *A_loc =
        (this->values.size() > 0) ? this->values.data() : nullptr;
      const NumberType *B_loc =
        (B.values.size() > 0) ? B.values.data() : nullptr;

      pgeadd(&trans_b,
             &n_rows,
             &n_columns,
             &beta,
             B_loc,
             &B.submatrix_row,
             &B.submatrix_column,
             B.descriptor,
             &alpha,
             A_loc,
             &submatrix_row,
             &submatrix_column,
             descriptor);
    }
  state = LAPACKSupport::matrix;
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::add(const NumberType                   a,
                                 const ScaLAPACKMatrix<NumberType> &B)
{
  add(B, 1, a, false);
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::Tadd(const NumberType                   a,
                                  const ScaLAPACKMatrix<NumberType> &B)
{
  add(B, 1, a, true);
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::mult(const NumberType                   b,
                                  const ScaLAPACKMatrix<NumberType> &B,
                                  const NumberType                   c,
                                  ScaLAPACKMatrix<NumberType> &      C,
                                  const bool transpose_A,
                                  const bool transpose_B) const
{
  Assert(this->grid == B.grid,
         ExcMessage("The matrices A and B need to have the same process grid"));
  Assert(C.grid == B.grid,
         ExcMessage("The matrices B and C need to have the same process grid"));

  // see for further info:
  // https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lgemm.htm
  if (!transpose_A && !transpose_B)
    {
      Assert(this->n_columns == B.n_rows,
             ExcDimensionMismatch(this->n_columns, B.n_rows));
      Assert(this->n_rows == C.n_rows,
             ExcDimensionMismatch(this->n_rows, C.n_rows));
      Assert(B.n_columns == C.n_columns,
             ExcDimensionMismatch(B.n_columns, C.n_columns));
      Assert(this->row_block_size == C.row_block_size,
             ExcDimensionMismatch(this->row_block_size, C.row_block_size));
      Assert(this->column_block_size == B.row_block_size,
             ExcDimensionMismatch(this->column_block_size, B.row_block_size));
      Assert(B.column_block_size == C.column_block_size,
             ExcDimensionMismatch(B.column_block_size, C.column_block_size));
    }
  else if (transpose_A && !transpose_B)
    {
      Assert(this->n_rows == B.n_rows,
             ExcDimensionMismatch(this->n_rows, B.n_rows));
      Assert(this->n_columns == C.n_rows,
             ExcDimensionMismatch(this->n_columns, C.n_rows));
      Assert(B.n_columns == C.n_columns,
             ExcDimensionMismatch(B.n_columns, C.n_columns));
      Assert(this->column_block_size == C.row_block_size,
             ExcDimensionMismatch(this->column_block_size, C.row_block_size));
      Assert(this->row_block_size == B.row_block_size,
             ExcDimensionMismatch(this->row_block_size, B.row_block_size));
      Assert(B.column_block_size == C.column_block_size,
             ExcDimensionMismatch(B.column_block_size, C.column_block_size));
    }
  else if (!transpose_A && transpose_B)
    {
      Assert(this->n_columns == B.n_columns,
             ExcDimensionMismatch(this->n_columns, B.n_columns));
      Assert(this->n_rows == C.n_rows,
             ExcDimensionMismatch(this->n_rows, C.n_rows));
      Assert(B.n_rows == C.n_columns,
             ExcDimensionMismatch(B.n_rows, C.n_columns));
      Assert(this->row_block_size == C.row_block_size,
             ExcDimensionMismatch(this->row_block_size, C.row_block_size));
      Assert(this->column_block_size == B.column_block_size,
             ExcDimensionMismatch(this->column_block_size,
                                  B.column_block_size));
      Assert(B.row_block_size == C.column_block_size,
             ExcDimensionMismatch(B.row_block_size, C.column_block_size));
    }
  else // if (transpose_A && transpose_B)
    {
      Assert(this->n_rows == B.n_columns,
             ExcDimensionMismatch(this->n_rows, B.n_columns));
      Assert(this->n_columns == C.n_rows,
             ExcDimensionMismatch(this->n_columns, C.n_rows));
      Assert(B.n_rows == C.n_columns,
             ExcDimensionMismatch(B.n_rows, C.n_columns));
      Assert(this->column_block_size == C.row_block_size,
             ExcDimensionMismatch(this->row_block_size, C.row_block_size));
      Assert(this->row_block_size == B.column_block_size,
             ExcDimensionMismatch(this->column_block_size, B.row_block_size));
      Assert(B.row_block_size == C.column_block_size,
             ExcDimensionMismatch(B.column_block_size, C.column_block_size));
    }

  if (this->grid->mpi_process_is_active)
    {
      char trans_a = transpose_A ? 'T' : 'N';
      char trans_b = transpose_B ? 'T' : 'N';

      const NumberType *A_loc =
        (this->values.size() > 0) ? this->values.data() : nullptr;
      const NumberType *B_loc =
        (B.values.size() > 0) ? B.values.data() : nullptr;
      NumberType *C_loc = (C.values.size() > 0) ? C.values.data() : nullptr;
      int         m     = C.n_rows;
      int         n     = C.n_columns;
      int         k     = transpose_A ? this->n_rows : this->n_columns;

      pgemm(&trans_a,
            &trans_b,
            &m,
            &n,
            &k,
            &b,
            A_loc,
            &(this->submatrix_row),
            &(this->submatrix_column),
            this->descriptor,
            B_loc,
            &B.submatrix_row,
            &B.submatrix_column,
            B.descriptor,
            &c,
            C_loc,
            &C.submatrix_row,
            &C.submatrix_column,
            C.descriptor);
    }
  C.state = LAPACKSupport::matrix;
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::mmult(ScaLAPACKMatrix<NumberType> &      C,
                                   const ScaLAPACKMatrix<NumberType> &B,
                                   const bool adding) const
{
  if (adding)
    mult(1., B, 1., C, false, false);
  else
    mult(1., B, 0, C, false, false);
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::Tmmult(ScaLAPACKMatrix<NumberType> &      C,
                                    const ScaLAPACKMatrix<NumberType> &B,
                                    const bool adding) const
{
  if (adding)
    mult(1., B, 1., C, true, false);
  else
    mult(1., B, 0, C, true, false);
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::mTmult(ScaLAPACKMatrix<NumberType> &      C,
                                    const ScaLAPACKMatrix<NumberType> &B,
                                    const bool adding) const
{
  if (adding)
    mult(1., B, 1., C, false, true);
  else
    mult(1., B, 0, C, false, true);
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::TmTmult(ScaLAPACKMatrix<NumberType> &      C,
                                     const ScaLAPACKMatrix<NumberType> &B,
                                     const bool adding) const
{
  if (adding)
    mult(1., B, 1., C, true, true);
  else
    mult(1., B, 0, C, true, true);
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::compute_cholesky_factorization()
{
  Assert(
    n_columns == n_rows && property == LAPACKSupport::Property::symmetric,
    ExcMessage(
      "Cholesky factorization can be applied to symmetric matrices only."));
  Assert(state == LAPACKSupport::matrix,
         ExcMessage(
           "Matrix has to be in Matrix state before calling this function."));

  if (grid->mpi_process_is_active)
    {
      int         info  = 0;
      NumberType *A_loc = this->values.data();
      // pdpotrf_(&uplo,&n_columns,A_loc,&submatrix_row,&submatrix_column,descriptor,&info);
      ppotrf(&uplo,
             &n_columns,
             A_loc,
             &submatrix_row,
             &submatrix_column,
             descriptor,
             &info);
      AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("ppotrf", info));
    }
  state    = LAPACKSupport::cholesky;
  property = (uplo == 'L' ? LAPACKSupport::lower_triangular :
                            LAPACKSupport::upper_triangular);
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::compute_lu_factorization()
{
  Assert(state == LAPACKSupport::matrix,
         ExcMessage(
           "Matrix has to be in Matrix state before calling this function."));

  if (grid->mpi_process_is_active)
    {
      int         info  = 0;
      NumberType *A_loc = this->values.data();

      const int iarow = indxg2p_(&submatrix_row,
                                 &row_block_size,
                                 &(grid->this_process_row),
                                 &first_process_row,
                                 &(grid->n_process_rows));
      const int mp    = numroc_(&n_rows,
                             &row_block_size,
                             &(grid->this_process_row),
                             &iarow,
                             &(grid->n_process_rows));
      ipiv.resize(mp + row_block_size);

      pgetrf(&n_rows,
             &n_columns,
             A_loc,
             &submatrix_row,
             &submatrix_column,
             descriptor,
             ipiv.data(),
             &info);
      AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("pgetrf", info));
    }
  state    = LAPACKSupport::State::lu;
  property = LAPACKSupport::Property::general;
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::invert()
{
  // Check whether matrix is symmetric and save flag.
  // If a Cholesky factorization has been applied previously,
  // the original matrix was symmetric.
  const bool is_symmetric = (property == LAPACKSupport::symmetric ||
                             state == LAPACKSupport::State::cholesky);

  // Check whether matrix is triangular and is in an unfactorized state.
  const bool is_triangular = (property == LAPACKSupport::upper_triangular ||
                              property == LAPACKSupport::lower_triangular) &&
                             (state == LAPACKSupport::State::matrix ||
                              state == LAPACKSupport::State::inverse_matrix);

  if (is_triangular)
    {
      if (grid->mpi_process_is_active)
        {
          const char uploTriangular =
            property == LAPACKSupport::upper_triangular ? 'U' : 'L';
          const char  diag  = 'N';
          int         info  = 0;
          NumberType *A_loc = this->values.data();
          ptrtri(&uploTriangular,
                 &diag,
                 &n_columns,
                 A_loc,
                 &submatrix_row,
                 &submatrix_column,
                 descriptor,
                 &info);
          AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("ptrtri", info));
          // The inversion is stored in the same part as the triangular matrix,
          // so we don't need to re-set the property here.
        }
    }
  else
    {
      // Matrix is neither in Cholesky nor LU state.
      // Compute the required factorizations based on the property of the
      // matrix.
      if (!(state == LAPACKSupport::State::lu ||
            state == LAPACKSupport::State::cholesky))
        {
          if (is_symmetric)
            compute_cholesky_factorization();
          else
            compute_lu_factorization();
        }
      if (grid->mpi_process_is_active)
        {
          int         info  = 0;
          NumberType *A_loc = this->values.data();

          if (is_symmetric)
            {
              ppotri(&uplo,
                     &n_columns,
                     A_loc,
                     &submatrix_row,
                     &submatrix_column,
                     descriptor,
                     &info);
              AssertThrow(info == 0,
                          LAPACKSupport::ExcErrorCode("ppotri", info));
              property = LAPACKSupport::Property::symmetric;
            }
          else
            {
              int lwork = -1, liwork = -1;
              work.resize(1);
              iwork.resize(1);

              pgetri(&n_columns,
                     A_loc,
                     &submatrix_row,
                     &submatrix_column,
                     descriptor,
                     ipiv.data(),
                     work.data(),
                     &lwork,
                     iwork.data(),
                     &liwork,
                     &info);

              AssertThrow(info == 0,
                          LAPACKSupport::ExcErrorCode("pgetri", info));
              lwork  = static_cast<int>(work[0]);
              liwork = iwork[0];
              work.resize(lwork);
              iwork.resize(liwork);

              pgetri(&n_columns,
                     A_loc,
                     &submatrix_row,
                     &submatrix_column,
                     descriptor,
                     ipiv.data(),
                     work.data(),
                     &lwork,
                     iwork.data(),
                     &liwork,
                     &info);

              AssertThrow(info == 0,
                          LAPACKSupport::ExcErrorCode("pgetri", info));
            }
        }
    }
  state = LAPACKSupport::State::inverse_matrix;
}



template <typename NumberType>
std::vector<NumberType>
ScaLAPACKMatrix<NumberType>::eigenpairs_symmetric_by_index(
  const std::pair<unsigned int, unsigned int> &index_limits,
  const bool                                   compute_eigenvectors)
{
  // check validity of index limits
  AssertIndexRange(index_limits.first, n_rows);
  AssertIndexRange(index_limits.second, n_rows);

  std::pair<unsigned int, unsigned int> idx =
    std::make_pair(std::min(index_limits.first, index_limits.second),
                   std::max(index_limits.first, index_limits.second));

  // compute all eigenvalues/eigenvectors
  if (idx.first == 0 && idx.second == static_cast<unsigned int>(n_rows - 1))
    return eigenpairs_symmetric(compute_eigenvectors);
  else
    return eigenpairs_symmetric(compute_eigenvectors, idx);
}



template <typename NumberType>
std::vector<NumberType>
ScaLAPACKMatrix<NumberType>::eigenpairs_symmetric_by_value(
  const std::pair<NumberType, NumberType> &value_limits,
  const bool                               compute_eigenvectors)
{
  Assert(!std::isnan(value_limits.first),
         ExcMessage("value_limits.first is NaN"));
  Assert(!std::isnan(value_limits.second),
         ExcMessage("value_limits.second is NaN"));

  std::pair<unsigned int, unsigned int> indices =
    std::make_pair(numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int);

  return eigenpairs_symmetric(compute_eigenvectors, indices, value_limits);
}



template <typename NumberType>
std::vector<NumberType>
ScaLAPACKMatrix<NumberType>::eigenpairs_symmetric(
  const bool                                   compute_eigenvectors,
  const std::pair<unsigned int, unsigned int> &eigenvalue_idx,
  const std::pair<NumberType, NumberType> &    eigenvalue_limits)
{
  Assert(state == LAPACKSupport::matrix,
         ExcMessage(
           "Matrix has to be in Matrix state before calling this function."));
  Assert(property == LAPACKSupport::symmetric,
         ExcMessage("Matrix has to be symmetric for this operation."));

  std::lock_guard<std::mutex> lock(mutex);

  const bool use_values = (std::isnan(eigenvalue_limits.first) ||
                           std::isnan(eigenvalue_limits.second)) ?
                            false :
                            true;
  const bool use_indices =
    ((eigenvalue_idx.first == numbers::invalid_unsigned_int) ||
     (eigenvalue_idx.second == numbers::invalid_unsigned_int)) ?
      false :
      true;

  Assert(
    !(use_values && use_indices),
    ExcMessage(
      "Prescribing both the index and value range for the eigenvalues is ambiguous"));

  // if computation of eigenvectors is not required use a sufficiently small
  // distributed matrix
  std::unique_ptr<ScaLAPACKMatrix<NumberType>> eigenvectors =
    compute_eigenvectors ?
      std::make_unique<ScaLAPACKMatrix<NumberType>>(n_rows,
                                                    grid,
                                                    row_block_size) :
      std::make_unique<ScaLAPACKMatrix<NumberType>>(
        grid->n_process_rows, grid->n_process_columns, grid, 1, 1);

  eigenvectors->property = property;
  // number of eigenvalues to be returned from psyevx; upon successful exit ev
  // contains the m seclected eigenvalues in ascending order set to all
  // eigenvaleus in case we will be using psyev.
  int                     m = n_rows;
  std::vector<NumberType> ev(n_rows);

  if (grid->mpi_process_is_active)
    {
      int info = 0;
      /*
       * for jobz==N only eigenvalues are computed, for jobz='V' also the
       * eigenvectors of the matrix are computed
       */
      char jobz  = compute_eigenvectors ? 'V' : 'N';
      char range = 'A';
      // default value is to compute all eigenvalues and optionally eigenvectors
      bool       all_eigenpairs = true;
      NumberType vl = NumberType(), vu = NumberType();
      int        il = 1, iu = 1;
      // number of eigenvectors to be returned;
      // upon successful exit the first m=nz columns contain the selected
      // eigenvectors (only if jobz=='V')
      int        nz     = 0;
      NumberType abstol = NumberType();

      // orfac decides which eigenvectors should be reorthogonalized
      // see
      // http://www.netlib.org/scalapack/explore-html/df/d1a/pdsyevx_8f_source.html
      // for explanation to keeps simple no reorthogonalized will be done by
      // setting orfac to 0
      NumberType orfac = 0;
      // contains the indices of eigenvectors that failed to converge
      std::vector<int> ifail;
      // This array contains indices of eigenvectors corresponding to
      // a cluster of eigenvalues that could not be reorthogonalized
      // due to insufficient workspace
      // see
      // http://www.netlib.org/scalapack/explore-html/df/d1a/pdsyevx_8f_source.html
      // for explanation
      std::vector<int> iclustr;
      // This array contains the gap between eigenvalues whose
      // eigenvectors could not be reorthogonalized.
      // see
      // http://www.netlib.org/scalapack/explore-html/df/d1a/pdsyevx_8f_source.html
      // for explanation
      std::vector<NumberType> gap(n_local_rows * n_local_columns);

      // index range for eigenvalues is not specified
      if (!use_indices)
        {
          // interval for eigenvalues is not specified and consequently all
          // eigenvalues/eigenpairs will be computed
          if (!use_values)
            {
              range          = 'A';
              all_eigenpairs = true;
            }
          else
            {
              range          = 'V';
              all_eigenpairs = false;
              vl = std::min(eigenvalue_limits.first, eigenvalue_limits.second);
              vu = std::max(eigenvalue_limits.first, eigenvalue_limits.second);
            }
        }
      else
        {
          range          = 'I';
          all_eigenpairs = false;
          // as Fortran starts counting/indexing from 1 unlike C/C++, where it
          // starts from 0
          il = std::min(eigenvalue_idx.first, eigenvalue_idx.second) + 1;
          iu = std::max(eigenvalue_idx.first, eigenvalue_idx.second) + 1;
        }
      NumberType *A_loc = this->values.data();
      /*
       * by setting lwork to -1 a workspace query for optimal length of work is
       * performed
       */
      int         lwork  = -1;
      int         liwork = -1;
      NumberType *eigenvectors_loc =
        (compute_eigenvectors ? eigenvectors->values.data() : nullptr);
      work.resize(1);
      iwork.resize(1);

      if (all_eigenpairs)
        {
          psyev(&jobz,
                &uplo,
                &n_rows,
                A_loc,
                &submatrix_row,
                &submatrix_column,
                descriptor,
                ev.data(),
                eigenvectors_loc,
                &eigenvectors->submatrix_row,
                &eigenvectors->submatrix_column,
                eigenvectors->descriptor,
                work.data(),
                &lwork,
                &info);
          AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("psyev", info));
        }
      else
        {
          char cmach = compute_eigenvectors ? 'U' : 'S';
          plamch(&(this->grid->blacs_context), &cmach, abstol);
          abstol *= 2;
          ifail.resize(n_rows);
          iclustr.resize(2 * grid->n_process_rows * grid->n_process_columns);
          gap.resize(grid->n_process_rows * grid->n_process_columns);

          psyevx(&jobz,
                 &range,
                 &uplo,
                 &n_rows,
                 A_loc,
                 &submatrix_row,
                 &submatrix_column,
                 descriptor,
                 &vl,
                 &vu,
                 &il,
                 &iu,
                 &abstol,
                 &m,
                 &nz,
                 ev.data(),
                 &orfac,
                 eigenvectors_loc,
                 &eigenvectors->submatrix_row,
                 &eigenvectors->submatrix_column,
                 eigenvectors->descriptor,
                 work.data(),
                 &lwork,
                 iwork.data(),
                 &liwork,
                 ifail.data(),
                 iclustr.data(),
                 gap.data(),
                 &info);
          AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("psyevx", info));
        }
      lwork = static_cast<int>(work[0]);
      work.resize(lwork);

      if (all_eigenpairs)
        {
          psyev(&jobz,
                &uplo,
                &n_rows,
                A_loc,
                &submatrix_row,
                &submatrix_column,
                descriptor,
                ev.data(),
                eigenvectors_loc,
                &eigenvectors->submatrix_row,
                &eigenvectors->submatrix_column,
                eigenvectors->descriptor,
                work.data(),
                &lwork,
                &info);

          AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("psyev", info));
        }
      else
        {
          liwork = iwork[0];
          AssertThrow(liwork > 0, ExcInternalError());
          iwork.resize(liwork);

          psyevx(&jobz,
                 &range,
                 &uplo,
                 &n_rows,
                 A_loc,
                 &submatrix_row,
                 &submatrix_column,
                 descriptor,
                 &vl,
                 &vu,
                 &il,
                 &iu,
                 &abstol,
                 &m,
                 &nz,
                 ev.data(),
                 &orfac,
                 eigenvectors_loc,
                 &eigenvectors->submatrix_row,
                 &eigenvectors->submatrix_column,
                 eigenvectors->descriptor,
                 work.data(),
                 &lwork,
                 iwork.data(),
                 &liwork,
                 ifail.data(),
                 iclustr.data(),
                 gap.data(),
                 &info);

          AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("psyevx", info));
        }
      // if eigenvectors are queried copy eigenvectors to original matrix
      // as the temporary matrix eigenvectors has identical dimensions and
      // block-cyclic distribution we simply swap the local array
      if (compute_eigenvectors)
        this->values.swap(eigenvectors->values);

      // adapt the size of ev to fit m upon return
      while (ev.size() > static_cast<size_type>(m))
        ev.pop_back();
    }
  /*
   * send number of computed eigenvalues to inactive processes
   */
  grid->send_to_inactive(&m, 1);

  /*
   * inactive processes have to resize array of eigenvalues
   */
  if (!grid->mpi_process_is_active)
    ev.resize(m);
  /*
   * send the eigenvalues to processors not being part of the process grid
   */
  grid->send_to_inactive(ev.data(), ev.size());

  /*
   * if only eigenvalues are queried the content of the matrix will be destroyed
   * if the eigenpairs are queried matrix A on exit stores the eigenvectors in
   * the columns
   */
  if (compute_eigenvectors)
    {
      property = LAPACKSupport::Property::general;
      state    = LAPACKSupport::eigenvalues;
    }
  else
    state = LAPACKSupport::unusable;

  return ev;
}



template <typename NumberType>
std::vector<NumberType>
ScaLAPACKMatrix<NumberType>::eigenpairs_symmetric_by_index_MRRR(
  const std::pair<unsigned int, unsigned int> &index_limits,
  const bool                                   compute_eigenvectors)
{
  // Check validity of index limits.
  AssertIndexRange(index_limits.first, static_cast<unsigned int>(n_rows));
  AssertIndexRange(index_limits.second, static_cast<unsigned int>(n_rows));

  const std::pair<unsigned int, unsigned int> idx =
    std::make_pair(std::min(index_limits.first, index_limits.second),
                   std::max(index_limits.first, index_limits.second));

  // Compute all eigenvalues/eigenvectors.
  if (idx.first == 0 && idx.second == static_cast<unsigned int>(n_rows - 1))
    return eigenpairs_symmetric_MRRR(compute_eigenvectors);
  else
    return eigenpairs_symmetric_MRRR(compute_eigenvectors, idx);
}



template <typename NumberType>
std::vector<NumberType>
ScaLAPACKMatrix<NumberType>::eigenpairs_symmetric_by_value_MRRR(
  const std::pair<NumberType, NumberType> &value_limits,
  const bool                               compute_eigenvectors)
{
  AssertIsFinite(value_limits.first);
  AssertIsFinite(value_limits.second);

  const std::pair<unsigned int, unsigned int> indices =
    std::make_pair(numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int);

  return eigenpairs_symmetric_MRRR(compute_eigenvectors, indices, value_limits);
}



template <typename NumberType>
std::vector<NumberType>
ScaLAPACKMatrix<NumberType>::eigenpairs_symmetric_MRRR(
  const bool                                   compute_eigenvectors,
  const std::pair<unsigned int, unsigned int> &eigenvalue_idx,
  const std::pair<NumberType, NumberType> &    eigenvalue_limits)
{
  Assert(state == LAPACKSupport::matrix,
         ExcMessage(
           "Matrix has to be in Matrix state before calling this function."));
  Assert(property == LAPACKSupport::symmetric,
         ExcMessage("Matrix has to be symmetric for this operation."));

  std::lock_guard<std::mutex> lock(mutex);

  const bool use_values = (std::isnan(eigenvalue_limits.first) ||
                           std::isnan(eigenvalue_limits.second)) ?
                            false :
                            true;
  const bool use_indices =
    ((eigenvalue_idx.first == numbers::invalid_unsigned_int) ||
     (eigenvalue_idx.second == numbers::invalid_unsigned_int)) ?
      false :
      true;

  Assert(
    !(use_values && use_indices),
    ExcMessage(
      "Prescribing both the index and value range for the eigenvalues is ambiguous"));

  // If computation of eigenvectors is not required, use a sufficiently small
  // distributed matrix.
  std::unique_ptr<ScaLAPACKMatrix<NumberType>> eigenvectors =
    compute_eigenvectors ?
      std::make_unique<ScaLAPACKMatrix<NumberType>>(n_rows,
                                                    grid,
                                                    row_block_size) :
      std::make_unique<ScaLAPACKMatrix<NumberType>>(
        grid->n_process_rows, grid->n_process_columns, grid, 1, 1);

  eigenvectors->property = property;
  // Number of eigenvalues to be returned from psyevr; upon successful exit ev
  // contains the m seclected eigenvalues in ascending order.
  int                     m = n_rows;
  std::vector<NumberType> ev(n_rows);

  // Number of eigenvectors to be returned;
  // Upon successful exit the first m=nz columns contain the selected
  // eigenvectors (only if jobz=='V').
  int nz = 0;

  if (grid->mpi_process_is_active)
    {
      int info = 0;
      /*
       * For jobz==N only eigenvalues are computed, for jobz='V' also the
       * eigenvectors of the matrix are computed.
       */
      char jobz = compute_eigenvectors ? 'V' : 'N';
      // Default value is to compute all eigenvalues and optionally
      // eigenvectors.
      char       range = 'A';
      NumberType vl = NumberType(), vu = NumberType();
      int        il = 1, iu = 1;

      // Index range for eigenvalues is not specified.
      if (!use_indices)
        {
          // Interval for eigenvalues is not specified and consequently all
          // eigenvalues/eigenpairs will be computed.
          if (!use_values)
            {
              range = 'A';
            }
          else
            {
              range = 'V';
              vl = std::min(eigenvalue_limits.first, eigenvalue_limits.second);
              vu = std::max(eigenvalue_limits.first, eigenvalue_limits.second);
            }
        }
      else
        {
          range = 'I';
          // As Fortran starts counting/indexing from 1 unlike C/C++, where it
          // starts from 0.
          il = std::min(eigenvalue_idx.first, eigenvalue_idx.second) + 1;
          iu = std::max(eigenvalue_idx.first, eigenvalue_idx.second) + 1;
        }
      NumberType *A_loc = this->values.data();

      /*
       * By setting lwork to -1 a workspace query for optimal length of work is
       * performed.
       */
      int         lwork  = -1;
      int         liwork = -1;
      NumberType *eigenvectors_loc =
        (compute_eigenvectors ? eigenvectors->values.data() : nullptr);
      work.resize(1);
      iwork.resize(1);

      psyevr(&jobz,
             &range,
             &uplo,
             &n_rows,
             A_loc,
             &submatrix_row,
             &submatrix_column,
             descriptor,
             &vl,
             &vu,
             &il,
             &iu,
             &m,
             &nz,
             ev.data(),
             eigenvectors_loc,
             &eigenvectors->submatrix_row,
             &eigenvectors->submatrix_column,
             eigenvectors->descriptor,
             work.data(),
             &lwork,
             iwork.data(),
             &liwork,
             &info);

      AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("psyevr", info));

      lwork = static_cast<int>(work[0]);
      work.resize(lwork);
      liwork = iwork[0];
      iwork.resize(liwork);

      psyevr(&jobz,
             &range,
             &uplo,
             &n_rows,
             A_loc,
             &submatrix_row,
             &submatrix_column,
             descriptor,
             &vl,
             &vu,
             &il,
             &iu,
             &m,
             &nz,
             ev.data(),
             eigenvectors_loc,
             &eigenvectors->submatrix_row,
             &eigenvectors->submatrix_column,
             eigenvectors->descriptor,
             work.data(),
             &lwork,
             iwork.data(),
             &liwork,
             &info);

      AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("psyevr", info));

      if (compute_eigenvectors)
        AssertThrow(
          m == nz,
          ExcMessage(
            "psyevr failed to compute all eigenvectors for the selected eigenvalues"));

      // If eigenvectors are queried, copy eigenvectors to original matrix.
      // As the temporary matrix eigenvectors has identical dimensions and
      // block-cyclic distribution we simply swap the local array.
      if (compute_eigenvectors)
        this->values.swap(eigenvectors->values);

      // Adapt the size of ev to fit m upon return.
      while (ev.size() > static_cast<size_type>(m))
        ev.pop_back();
    }
  /*
   * Send number of computed eigenvalues to inactive processes.
   */
  grid->send_to_inactive(&m, 1);

  /*
   * Inactive processes have to resize array of eigenvalues.
   */
  if (!grid->mpi_process_is_active)
    ev.resize(m);
  /*
   * Send the eigenvalues to processors not being part of the process grid.
   */
  grid->send_to_inactive(ev.data(), ev.size());

  /*
   * If only eigenvalues are queried, the content of the matrix will be
   * destroyed. If the eigenpairs are queried, matrix A on exit stores the
   * eigenvectors in the columns.
   */
  if (compute_eigenvectors)
    {
      property = LAPACKSupport::Property::general;
      state    = LAPACKSupport::eigenvalues;
    }
  else
    state = LAPACKSupport::unusable;

  return ev;
}



template <typename NumberType>
std::vector<NumberType>
ScaLAPACKMatrix<NumberType>::compute_SVD(ScaLAPACKMatrix<NumberType> *U,
                                         ScaLAPACKMatrix<NumberType> *VT)
{
  Assert(state == LAPACKSupport::matrix,
         ExcMessage(
           "Matrix has to be in Matrix state before calling this function."));
  Assert(row_block_size == column_block_size,
         ExcDimensionMismatch(row_block_size, column_block_size));

  const bool left_singluar_vectors  = (U != nullptr) ? true : false;
  const bool right_singluar_vectors = (VT != nullptr) ? true : false;

  if (left_singluar_vectors)
    {
      Assert(n_rows == U->n_rows, ExcDimensionMismatch(n_rows, U->n_rows));
      Assert(U->n_rows == U->n_columns,
             ExcDimensionMismatch(U->n_rows, U->n_columns));
      Assert(row_block_size == U->row_block_size,
             ExcDimensionMismatch(row_block_size, U->row_block_size));
      Assert(column_block_size == U->column_block_size,
             ExcDimensionMismatch(column_block_size, U->column_block_size));
      Assert(grid->blacs_context == U->grid->blacs_context,
             ExcDimensionMismatch(grid->blacs_context, U->grid->blacs_context));
    }
  if (right_singluar_vectors)
    {
      Assert(n_columns == VT->n_rows,
             ExcDimensionMismatch(n_columns, VT->n_rows));
      Assert(VT->n_rows == VT->n_columns,
             ExcDimensionMismatch(VT->n_rows, VT->n_columns));
      Assert(row_block_size == VT->row_block_size,
             ExcDimensionMismatch(row_block_size, VT->row_block_size));
      Assert(column_block_size == VT->column_block_size,
             ExcDimensionMismatch(column_block_size, VT->column_block_size));
      Assert(grid->blacs_context == VT->grid->blacs_context,
             ExcDimensionMismatch(grid->blacs_context,
                                  VT->grid->blacs_context));
    }
  std::lock_guard<std::mutex> lock(mutex);

  std::vector<NumberType> sv(std::min(n_rows, n_columns));

  if (grid->mpi_process_is_active)
    {
      char        jobu   = left_singluar_vectors ? 'V' : 'N';
      char        jobvt  = right_singluar_vectors ? 'V' : 'N';
      NumberType *A_loc  = this->values.data();
      NumberType *U_loc  = left_singluar_vectors ? U->values.data() : nullptr;
      NumberType *VT_loc = right_singluar_vectors ? VT->values.data() : nullptr;
      int         info   = 0;
      /*
       * by setting lwork to -1 a workspace query for optimal length of work is
       * performed
       */
      int lwork = -1;
      work.resize(1);

      pgesvd(&jobu,
             &jobvt,
             &n_rows,
             &n_columns,
             A_loc,
             &submatrix_row,
             &submatrix_column,
             descriptor,
             &*sv.begin(),
             U_loc,
             &U->submatrix_row,
             &U->submatrix_column,
             U->descriptor,
             VT_loc,
             &VT->submatrix_row,
             &VT->submatrix_column,
             VT->descriptor,
             work.data(),
             &lwork,
             &info);
      AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("pgesvd", info));

      lwork = static_cast<int>(work[0]);
      work.resize(lwork);

      pgesvd(&jobu,
             &jobvt,
             &n_rows,
             &n_columns,
             A_loc,
             &submatrix_row,
             &submatrix_column,
             descriptor,
             &*sv.begin(),
             U_loc,
             &U->submatrix_row,
             &U->submatrix_column,
             U->descriptor,
             VT_loc,
             &VT->submatrix_row,
             &VT->submatrix_column,
             VT->descriptor,
             work.data(),
             &lwork,
             &info);
      AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("pgesvd", info));
    }

  /*
   * send the singular values to processors not being part of the process grid
   */
  grid->send_to_inactive(sv.data(), sv.size());

  property = LAPACKSupport::Property::general;
  state    = LAPACKSupport::State::unusable;

  return sv;
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::least_squares(ScaLAPACKMatrix<NumberType> &B,
                                           const bool transpose)
{
  Assert(grid == B.grid,
         ExcMessage("The matrices A and B need to have the same process grid"));
  Assert(state == LAPACKSupport::matrix,
         ExcMessage(
           "Matrix has to be in Matrix state before calling this function."));
  Assert(B.state == LAPACKSupport::matrix,
         ExcMessage(
           "Matrix B has to be in Matrix state before calling this function."));

  if (transpose)
    {
      Assert(n_columns == B.n_rows, ExcDimensionMismatch(n_columns, B.n_rows));
    }
  else
    {
      Assert(n_rows == B.n_rows, ExcDimensionMismatch(n_rows, B.n_rows));
    }

  // see
  // https://www.ibm.com/support/knowledgecenter/en/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lgels.htm
  Assert(row_block_size == column_block_size,
         ExcMessage(
           "Use identical block sizes for rows and columns of matrix A"));
  Assert(B.row_block_size == B.column_block_size,
         ExcMessage(
           "Use identical block sizes for rows and columns of matrix B"));
  Assert(row_block_size == B.row_block_size,
         ExcMessage(
           "Use identical block-cyclic distribution for matrices A and B"));

  std::lock_guard<std::mutex> lock(mutex);

  if (grid->mpi_process_is_active)
    {
      char        trans = transpose ? 'T' : 'N';
      NumberType *A_loc = this->values.data();
      NumberType *B_loc = B.values.data();
      int         info  = 0;
      /*
       * by setting lwork to -1 a workspace query for optimal length of work is
       * performed
       */
      int lwork = -1;
      work.resize(1);

      pgels(&trans,
            &n_rows,
            &n_columns,
            &B.n_columns,
            A_loc,
            &submatrix_row,
            &submatrix_column,
            descriptor,
            B_loc,
            &B.submatrix_row,
            &B.submatrix_column,
            B.descriptor,
            work.data(),
            &lwork,
            &info);
      AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("pgels", info));

      lwork = static_cast<int>(work[0]);
      work.resize(lwork);

      pgels(&trans,
            &n_rows,
            &n_columns,
            &B.n_columns,
            A_loc,
            &submatrix_row,
            &submatrix_column,
            descriptor,
            B_loc,
            &B.submatrix_row,
            &B.submatrix_column,
            B.descriptor,
            work.data(),
            &lwork,
            &info);
      AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("pgels", info));
    }
  state = LAPACKSupport::State::unusable;
}



template <typename NumberType>
unsigned int
ScaLAPACKMatrix<NumberType>::pseudoinverse(const NumberType ratio)
{
  Assert(state == LAPACKSupport::matrix,
         ExcMessage(
           "Matrix has to be in Matrix state before calling this function."));
  Assert(row_block_size == column_block_size,
         ExcMessage(
           "Use identical block sizes for rows and columns of matrix A"));
  Assert(
    ratio > 0. && ratio < 1.,
    ExcMessage(
      "input parameter ratio has to be larger than zero and smaller than 1"));

  ScaLAPACKMatrix<NumberType> U(n_rows,
                                n_rows,
                                grid,
                                row_block_size,
                                row_block_size,
                                LAPACKSupport::Property::general);
  ScaLAPACKMatrix<NumberType> VT(n_columns,
                                 n_columns,
                                 grid,
                                 row_block_size,
                                 row_block_size,
                                 LAPACKSupport::Property::general);
  std::vector<NumberType>     sv = this->compute_SVD(&U, &VT);
  AssertThrow(sv[0] > std::numeric_limits<NumberType>::min(),
              ExcMessage("Matrix has rank 0"));

  // Get number of singular values fulfilling the following: sv[i] > sv[0] *
  // ratio Obviously, 0-th element already satisfies sv[0] > sv[0] * ratio The
  // singular values in sv are ordered by descending value so we break out of
  // the loop if a singular value is smaller than sv[0] * ratio.
  unsigned int            n_sv = 1;
  std::vector<NumberType> inv_sigma;
  inv_sigma.push_back(1 / sv[0]);

  for (unsigned int i = 1; i < sv.size(); ++i)
    if (sv[i] > sv[0] * ratio)
      {
        ++n_sv;
        inv_sigma.push_back(1 / sv[i]);
      }
    else
      break;

  // For the matrix multiplication we use only the columns of U and rows of VT
  // which are associated with singular values larger than the limit. That saves
  // computational time for matrices with rank significantly smaller than
  // min(n_rows,n_columns)
  ScaLAPACKMatrix<NumberType> U_R(n_rows,
                                  n_sv,
                                  grid,
                                  row_block_size,
                                  row_block_size,
                                  LAPACKSupport::Property::general);
  ScaLAPACKMatrix<NumberType> VT_R(n_sv,
                                   n_columns,
                                   grid,
                                   row_block_size,
                                   row_block_size,
                                   LAPACKSupport::Property::general);
  U.copy_to(U_R,
            std::make_pair(0, 0),
            std::make_pair(0, 0),
            std::make_pair(n_rows, n_sv));
  VT.copy_to(VT_R,
             std::make_pair(0, 0),
             std::make_pair(0, 0),
             std::make_pair(n_sv, n_columns));

  VT_R.scale_rows(inv_sigma);
  this->reinit(n_columns,
               n_rows,
               this->grid,
               row_block_size,
               column_block_size,
               LAPACKSupport::Property::general);
  VT_R.mult(1, U_R, 0, *this, true, true);
  state = LAPACKSupport::State::inverse_matrix;
  return n_sv;
}



template <typename NumberType>
NumberType
ScaLAPACKMatrix<NumberType>::reciprocal_condition_number(
  const NumberType a_norm) const
{
  Assert(state == LAPACKSupport::cholesky,
         ExcMessage(
           "Matrix has to be in Cholesky state before calling this function."));
  std::lock_guard<std::mutex> lock(mutex);
  NumberType                  rcond = 0.;

  if (grid->mpi_process_is_active)
    {
      int liwork = n_local_rows;
      iwork.resize(liwork);

      int               info  = 0;
      const NumberType *A_loc = this->values.data();

      // by setting lwork to -1 a workspace query for optimal length of work is
      // performed
      int lwork = -1;
      work.resize(1);
      ppocon(&uplo,
             &n_columns,
             A_loc,
             &submatrix_row,
             &submatrix_column,
             descriptor,
             &a_norm,
             &rcond,
             work.data(),
             &lwork,
             iwork.data(),
             &liwork,
             &info);
      AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("pdpocon", info));
      lwork = static_cast<int>(std::ceil(work[0]));
      work.resize(lwork);

      // now the actual run:
      ppocon(&uplo,
             &n_columns,
             A_loc,
             &submatrix_row,
             &submatrix_column,
             descriptor,
             &a_norm,
             &rcond,
             work.data(),
             &lwork,
             iwork.data(),
             &liwork,
             &info);
      AssertThrow(info == 0, LAPACKSupport::ExcErrorCode("pdpocon", info));
    }
  grid->send_to_inactive(&rcond);
  return rcond;
}



template <typename NumberType>
NumberType
ScaLAPACKMatrix<NumberType>::l1_norm() const
{
  const char type('O');

  if (property == LAPACKSupport::symmetric)
    return norm_symmetric(type);
  else
    return norm_general(type);
}



template <typename NumberType>
NumberType
ScaLAPACKMatrix<NumberType>::linfty_norm() const
{
  const char type('I');

  if (property == LAPACKSupport::symmetric)
    return norm_symmetric(type);
  else
    return norm_general(type);
}



template <typename NumberType>
NumberType
ScaLAPACKMatrix<NumberType>::frobenius_norm() const
{
  const char type('F');

  if (property == LAPACKSupport::symmetric)
    return norm_symmetric(type);
  else
    return norm_general(type);
}



template <typename NumberType>
NumberType
ScaLAPACKMatrix<NumberType>::norm_general(const char type) const
{
  Assert(state == LAPACKSupport::matrix ||
           state == LAPACKSupport::inverse_matrix,
         ExcMessage("norms can be called in matrix state only."));
  std::lock_guard<std::mutex> lock(mutex);
  NumberType                  res = 0.;

  if (grid->mpi_process_is_active)
    {
      const int iarow = indxg2p_(&submatrix_row,
                                 &row_block_size,
                                 &(grid->this_process_row),
                                 &first_process_row,
                                 &(grid->n_process_rows));
      const int iacol = indxg2p_(&submatrix_column,
                                 &column_block_size,
                                 &(grid->this_process_column),
                                 &first_process_column,
                                 &(grid->n_process_columns));
      const int mp0   = numroc_(&n_rows,
                              &row_block_size,
                              &(grid->this_process_row),
                              &iarow,
                              &(grid->n_process_rows));
      const int nq0   = numroc_(&n_columns,
                              &column_block_size,
                              &(grid->this_process_column),
                              &iacol,
                              &(grid->n_process_columns));

      // type='M': compute largest absolute value
      // type='F' || type='E': compute Frobenius norm
      // type='0' || type='1': compute infinity norm
      int lwork = 0; // for type == 'M' || type == 'F' || type == 'E'
      if (type == 'O' || type == '1')
        lwork = nq0;
      else if (type == 'I')
        lwork = mp0;

      work.resize(lwork);
      const NumberType *A_loc = this->values.begin();
      res                     = plange(&type,
                   &n_rows,
                   &n_columns,
                   A_loc,
                   &submatrix_row,
                   &submatrix_column,
                   descriptor,
                   work.data());
    }
  grid->send_to_inactive(&res);
  return res;
}



template <typename NumberType>
NumberType
ScaLAPACKMatrix<NumberType>::norm_symmetric(const char type) const
{
  Assert(state == LAPACKSupport::matrix ||
           state == LAPACKSupport::inverse_matrix,
         ExcMessage("norms can be called in matrix state only."));
  Assert(property == LAPACKSupport::symmetric,
         ExcMessage("Matrix has to be symmetric for this operation."));
  std::lock_guard<std::mutex> lock(mutex);
  NumberType                  res = 0.;

  if (grid->mpi_process_is_active)
    {
      // int IROFFA = MOD( IA-1, MB_A )
      // int ICOFFA = MOD( JA-1, NB_A )
      const int lcm =
        ilcm_(&(grid->n_process_rows), &(grid->n_process_columns));
      const int v2 = lcm / (grid->n_process_rows);

      const int IAROW = indxg2p_(&submatrix_row,
                                 &row_block_size,
                                 &(grid->this_process_row),
                                 &first_process_row,
                                 &(grid->n_process_rows));
      const int IACOL = indxg2p_(&submatrix_column,
                                 &column_block_size,
                                 &(grid->this_process_column),
                                 &first_process_column,
                                 &(grid->n_process_columns));
      const int Np0   = numroc_(&n_columns /*+IROFFA*/,
                              &row_block_size,
                              &(grid->this_process_row),
                              &IAROW,
                              &(grid->n_process_rows));
      const int Nq0   = numroc_(&n_columns /*+ICOFFA*/,
                              &column_block_size,
                              &(grid->this_process_column),
                              &IACOL,
                              &(grid->n_process_columns));

      const int v1  = iceil_(&Np0, &row_block_size);
      const int ldw = (n_local_rows == n_local_columns) ?
                        0 :
                        row_block_size * iceil_(&v1, &v2);

      const int lwork =
        (type == 'M' || type == 'F' || type == 'E') ? 0 : 2 * Nq0 + Np0 + ldw;
      work.resize(lwork);
      const NumberType *A_loc = this->values.begin();
      res                     = plansy(&type,
                   &uplo,
                   &n_columns,
                   A_loc,
                   &submatrix_row,
                   &submatrix_column,
                   descriptor,
                   work.data());
    }
  grid->send_to_inactive(&res);
  return res;
}



#  ifdef DEAL_II_WITH_HDF5
namespace internal
{
  namespace
  {
    void
    create_HDF5_state_enum_id(hid_t &state_enum_id)
    {
      // create HDF5 enum type for LAPACKSupport::State
      LAPACKSupport::State val;
      state_enum_id = H5Tcreate(H5T_ENUM, sizeof(LAPACKSupport::State));
      val           = LAPACKSupport::State::cholesky;
      herr_t status = H5Tenum_insert(state_enum_id, "cholesky", &val);
      AssertThrow(status >= 0, ExcInternalError());
      val    = LAPACKSupport::State::eigenvalues;
      status = H5Tenum_insert(state_enum_id, "eigenvalues", &val);
      AssertThrow(status >= 0, ExcInternalError());
      val    = LAPACKSupport::State::inverse_matrix;
      status = H5Tenum_insert(state_enum_id, "inverse_matrix", &val);
      AssertThrow(status >= 0, ExcInternalError());
      val    = LAPACKSupport::State::inverse_svd;
      status = H5Tenum_insert(state_enum_id, "inverse_svd", &val);
      AssertThrow(status >= 0, ExcInternalError());
      val    = LAPACKSupport::State::lu;
      status = H5Tenum_insert(state_enum_id, "lu", &val);
      AssertThrow(status >= 0, ExcInternalError());
      val    = LAPACKSupport::State::matrix;
      status = H5Tenum_insert(state_enum_id, "matrix", &val);
      AssertThrow(status >= 0, ExcInternalError());
      val    = LAPACKSupport::State::svd;
      status = H5Tenum_insert(state_enum_id, "svd", &val);
      AssertThrow(status >= 0, ExcInternalError());
      val    = LAPACKSupport::State::unusable;
      status = H5Tenum_insert(state_enum_id, "unusable", &val);
      AssertThrow(status >= 0, ExcInternalError());
    }

    void
    create_HDF5_property_enum_id(hid_t &property_enum_id)
    {
      // create HDF5 enum type for LAPACKSupport::Property
      property_enum_id = H5Tcreate(H5T_ENUM, sizeof(LAPACKSupport::Property));
      LAPACKSupport::Property prop = LAPACKSupport::Property::diagonal;
      herr_t status = H5Tenum_insert(property_enum_id, "diagonal", &prop);
      AssertThrow(status >= 0, ExcInternalError());
      prop   = LAPACKSupport::Property::general;
      status = H5Tenum_insert(property_enum_id, "general", &prop);
      AssertThrow(status >= 0, ExcInternalError());
      prop   = LAPACKSupport::Property::hessenberg;
      status = H5Tenum_insert(property_enum_id, "hessenberg", &prop);
      AssertThrow(status >= 0, ExcInternalError());
      prop   = LAPACKSupport::Property::lower_triangular;
      status = H5Tenum_insert(property_enum_id, "lower_triangular", &prop);
      AssertThrow(status >= 0, ExcInternalError());
      prop   = LAPACKSupport::Property::symmetric;
      status = H5Tenum_insert(property_enum_id, "symmetric", &prop);
      AssertThrow(status >= 0, ExcInternalError());
      prop   = LAPACKSupport::Property::upper_triangular;
      status = H5Tenum_insert(property_enum_id, "upper_triangular", &prop);
      AssertThrow(status >= 0, ExcInternalError());
    }
  } // namespace
} // namespace internal
#  endif



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::save(
  const std::string &                          filename,
  const std::pair<unsigned int, unsigned int> &chunk_size) const
{
#  ifndef DEAL_II_WITH_HDF5
  (void)filename;
  (void)chunk_size;
  AssertThrow(false, ExcMessage("HDF5 support is disabled."));
#  else

  std::pair<unsigned int, unsigned int> chunks_size_ = chunk_size;

  if (chunks_size_.first == numbers::invalid_unsigned_int ||
      chunks_size_.second == numbers::invalid_unsigned_int)
    {
      // default: store the matrix in chunks of columns
      chunks_size_.first  = n_rows;
      chunks_size_.second = 1;
    }
  Assert(chunks_size_.first > 0,
         ExcMessage("The row chunk size must be larger than 0."));
  AssertIndexRange(chunks_size_.first, n_rows + 1);
  Assert(chunks_size_.second > 0,
         ExcMessage("The column chunk size must be larger than 0."));
  AssertIndexRange(chunks_size_.second, n_columns + 1);

#    ifdef H5_HAVE_PARALLEL
  // implementation for configurations equipped with a parallel file system
  save_parallel(filename, chunks_size_);

#    else
  // implementation for configurations with no parallel file system
  save_serial(filename, chunks_size_);

#    endif
#  endif
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::save_serial(
  const std::string &                          filename,
  const std::pair<unsigned int, unsigned int> &chunk_size) const
{
#  ifndef DEAL_II_WITH_HDF5
  (void)filename;
  (void)chunk_size;
  Assert(false, ExcInternalError());
#  else

  /*
   * The content of the distributed matrix is copied to a matrix using a 1x1
   * process grid. Therefore, one process has all the data and can write it to a
   * file.
   *
   * Create a 1x1 column grid which will be used to initialize
   * an effectively serial ScaLAPACK matrix to gather the contents from the
   * current object
   */
  const auto column_grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(this->grid->mpi_communicator,
                                                  1,
                                                  1);

  const int                   MB = n_rows, NB = n_columns;
  ScaLAPACKMatrix<NumberType> tmp(n_rows, n_columns, column_grid, MB, NB);
  copy_to(tmp);

  // the 1x1 grid has only one process and this one writes
  // the content of the matrix to the HDF5 file
  if (tmp.grid->mpi_process_is_active)
    {
      herr_t status;

      // create a new file using default properties
      hid_t file_id =
        H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      // modify dataset creation properties, i.e. enable chunking
      hsize_t chunk_dims[2];
      // revert order of rows and columns as ScaLAPACK uses column-major
      // ordering
      chunk_dims[0]       = chunk_size.second;
      chunk_dims[1]       = chunk_size.first;
      hid_t data_property = H5Pcreate(H5P_DATASET_CREATE);
      status              = H5Pset_chunk(data_property, 2, chunk_dims);
      AssertThrow(status >= 0, ExcIO());

      // create the data space for the dataset
      hsize_t dims[2];
      // change order of rows and columns as ScaLAPACKMatrix uses column major
      // ordering
      dims[0]            = n_columns;
      dims[1]            = n_rows;
      hid_t dataspace_id = H5Screate_simple(2, dims, nullptr);

      // create the dataset within the file using chunk creation properties
      hid_t type_id    = hdf5_type_id(tmp.values.data());
      hid_t dataset_id = H5Dcreate2(file_id,
                                    "/matrix",
                                    type_id,
                                    dataspace_id,
                                    H5P_DEFAULT,
                                    data_property,
                                    H5P_DEFAULT);

      // write the dataset
      status = H5Dwrite(
        dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp.values.data());
      AssertThrow(status >= 0, ExcIO());

      // create HDF5 enum type for LAPACKSupport::State and
      // LAPACKSupport::Property
      hid_t state_enum_id, property_enum_id;
      internal::create_HDF5_state_enum_id(state_enum_id);
      internal::create_HDF5_property_enum_id(property_enum_id);

      // create the data space for the state enum
      hsize_t dims_state[1];
      dims_state[0]              = 1;
      hid_t state_enum_dataspace = H5Screate_simple(1, dims_state, nullptr);
      // create the dataset for the state enum
      hid_t state_enum_dataset = H5Dcreate2(file_id,
                                            "/state",
                                            state_enum_id,
                                            state_enum_dataspace,
                                            H5P_DEFAULT,
                                            H5P_DEFAULT,
                                            H5P_DEFAULT);
      // write the dataset for the state enum
      status = H5Dwrite(state_enum_dataset,
                        state_enum_id,
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        &state);
      AssertThrow(status >= 0, ExcIO());

      // create the data space for the property enum
      hsize_t dims_property[1];
      dims_property[0] = 1;
      hid_t property_enum_dataspace =
        H5Screate_simple(1, dims_property, nullptr);
      // create the dataset for the property enum
      hid_t property_enum_dataset = H5Dcreate2(file_id,
                                               "/property",
                                               property_enum_id,
                                               property_enum_dataspace,
                                               H5P_DEFAULT,
                                               H5P_DEFAULT,
                                               H5P_DEFAULT);
      // write the dataset for the property enum
      status = H5Dwrite(property_enum_dataset,
                        property_enum_id,
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        &property);
      AssertThrow(status >= 0, ExcIO());

      // end access to the datasets and release resources used by them
      status = H5Dclose(dataset_id);
      AssertThrow(status >= 0, ExcIO());
      status = H5Dclose(state_enum_dataset);
      AssertThrow(status >= 0, ExcIO());
      status = H5Dclose(property_enum_dataset);
      AssertThrow(status >= 0, ExcIO());

      // terminate access to the data spaces
      status = H5Sclose(dataspace_id);
      AssertThrow(status >= 0, ExcIO());
      status = H5Sclose(state_enum_dataspace);
      AssertThrow(status >= 0, ExcIO());
      status = H5Sclose(property_enum_dataspace);
      AssertThrow(status >= 0, ExcIO());

      // release enum data types
      status = H5Tclose(state_enum_id);
      AssertThrow(status >= 0, ExcIO());
      status = H5Tclose(property_enum_id);
      AssertThrow(status >= 0, ExcIO());

      // release the creation property
      status = H5Pclose(data_property);
      AssertThrow(status >= 0, ExcIO());

      // close the file.
      status = H5Fclose(file_id);
      AssertThrow(status >= 0, ExcIO());
    }
#  endif
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::save_parallel(
  const std::string &                          filename,
  const std::pair<unsigned int, unsigned int> &chunk_size) const
{
#  ifndef DEAL_II_WITH_HDF5
  (void)filename;
  (void)chunk_size;
  Assert(false, ExcInternalError());
#  else

  const unsigned int n_mpi_processes(
    Utilities::MPI::n_mpi_processes(this->grid->mpi_communicator));
  MPI_Info info = MPI_INFO_NULL;
  /*
   * The content of the distributed matrix is copied to a matrix using a
   * 1xn_processes process grid. Therefore, the processes hold contiguous chunks
   * of the matrix, which they can write to the file
   *
   * Create a 1xn_processes column grid
   */
  const auto column_grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(this->grid->mpi_communicator,
                                                  1,
                                                  n_mpi_processes);

  const int MB = n_rows;
  /*
   * If the ratio n_columns/n_mpi_processes is smaller than the column block
   * size of the original matrix, the redistribution and saving of the matrix
   * requires a significant amount of MPI communication. Therefore, it is better
   * to set a minimum value for the block size NB, causing only
   * ceil(n_columns/NB) processes being actively involved in saving the matrix.
   * Example: A 2*10^9 x 400 matrix is distributed on a 80 x 5 process grid
   * using block size 32. Instead of distributing the matrix on a 1 x 400
   * process grid with a row block size of 2*10^9 and a column block size of 1,
   * the minimum value for NB yields that only ceil(400/32)=13 processes will be
   * writing the matrix to disk.
   */
  const int NB = std::max(static_cast<int>(std::ceil(
                            static_cast<double>(n_columns) / n_mpi_processes)),
                          column_block_size);

  ScaLAPACKMatrix<NumberType> tmp(n_rows, n_columns, column_grid, MB, NB);
  copy_to(tmp);

  // get pointer to data held by the process
  NumberType *data = (tmp.values.size() > 0) ? tmp.values.data() : nullptr;

  herr_t status;
  // dataset dimensions
  hsize_t dims[2];

  // set up file access property list with parallel I/O access
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  status         = H5Pset_fapl_mpio(plist_id, tmp.grid->mpi_communicator, info);
  AssertThrow(status >= 0, ExcIO());

  // create a new file collectively and release property list identifier
  hid_t file_id =
    H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  status = H5Pclose(plist_id);
  AssertThrow(status >= 0, ExcIO());

  // As ScaLAPACK, and therefore the class ScaLAPACKMatrix, uses column-major
  // ordering but HDF5 row-major ordering, we have to reverse entries related to
  // columns and rows in the following. create the dataspace for the dataset
  dims[0] = tmp.n_columns;
  dims[1] = tmp.n_rows;

  hid_t filespace = H5Screate_simple(2, dims, nullptr);

  // create the chunked dataset with default properties and close filespace
  hsize_t chunk_dims[2];
  // revert order of rows and columns as ScaLAPACK uses column-major ordering
  chunk_dims[0] = chunk_size.second;
  chunk_dims[1] = chunk_size.first;
  plist_id      = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(plist_id, 2, chunk_dims);
  hid_t type_id = hdf5_type_id(data);
  hid_t dset_id = H5Dcreate2(
    file_id, "/matrix", type_id, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);

  status = H5Sclose(filespace);
  AssertThrow(status >= 0, ExcIO());

  status = H5Pclose(plist_id);
  AssertThrow(status >= 0, ExcIO());

  // gather the number of local rows and columns from all processes
  std::vector<int> proc_n_local_rows(n_mpi_processes),
    proc_n_local_columns(n_mpi_processes);
  int ierr = MPI_Allgather(&tmp.n_local_rows,
                           1,
                           MPI_INT,
                           proc_n_local_rows.data(),
                           1,
                           MPI_INT,
                           tmp.grid->mpi_communicator);
  AssertThrowMPI(ierr);
  ierr = MPI_Allgather(&tmp.n_local_columns,
                       1,
                       MPI_INT,
                       proc_n_local_columns.data(),
                       1,
                       MPI_INT,
                       tmp.grid->mpi_communicator);
  AssertThrowMPI(ierr);

  const unsigned int my_rank(
    Utilities::MPI::this_mpi_process(tmp.grid->mpi_communicator));

  // hyperslab selection parameters
  // each process defines dataset in memory and writes it to the hyperslab in
  // the file
  hsize_t count[2];
  count[0]       = tmp.n_local_columns;
  count[1]       = tmp.n_rows;
  hid_t memspace = H5Screate_simple(2, count, nullptr);

  hsize_t offset[2] = {0};
  for (unsigned int i = 0; i < my_rank; ++i)
    offset[0] += proc_n_local_columns[i];

  // select hyperslab in the file.
  filespace = H5Dget_space(dset_id);
  status    = H5Sselect_hyperslab(
    filespace, H5S_SELECT_SET, offset, nullptr, count, nullptr);
  AssertThrow(status >= 0, ExcIO());

  // create property list for independent dataset write
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  status   = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  AssertThrow(status >= 0, ExcIO());

  // process with no data will not participate in writing to the file
  if (tmp.values.size() > 0)
    {
      status = H5Dwrite(dset_id, type_id, memspace, filespace, plist_id, data);
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

  // before writing the state and property to file wait for
  // all processes to finish writing the matrix content to the file
  ierr = MPI_Barrier(tmp.grid->mpi_communicator);
  AssertThrowMPI(ierr);

  // only root process will write state and property to the file
  if (tmp.grid->this_mpi_process == 0)
    {
      // open file using default properties
      hid_t file_id_reopen =
        H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

      // create HDF5 enum type for LAPACKSupport::State and
      // LAPACKSupport::Property
      hid_t state_enum_id, property_enum_id;
      internal::create_HDF5_state_enum_id(state_enum_id);
      internal::create_HDF5_property_enum_id(property_enum_id);

      // create the data space for the state enum
      hsize_t dims_state[1];
      dims_state[0]              = 1;
      hid_t state_enum_dataspace = H5Screate_simple(1, dims_state, nullptr);
      // create the dataset for the state enum
      hid_t state_enum_dataset = H5Dcreate2(file_id_reopen,
                                            "/state",
                                            state_enum_id,
                                            state_enum_dataspace,
                                            H5P_DEFAULT,
                                            H5P_DEFAULT,
                                            H5P_DEFAULT);
      // write the dataset for the state enum
      status = H5Dwrite(state_enum_dataset,
                        state_enum_id,
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        &state);
      AssertThrow(status >= 0, ExcIO());

      // create the data space for the property enum
      hsize_t dims_property[1];
      dims_property[0] = 1;
      hid_t property_enum_dataspace =
        H5Screate_simple(1, dims_property, nullptr);
      // create the dataset for the property enum
      hid_t property_enum_dataset = H5Dcreate2(file_id_reopen,
                                               "/property",
                                               property_enum_id,
                                               property_enum_dataspace,
                                               H5P_DEFAULT,
                                               H5P_DEFAULT,
                                               H5P_DEFAULT);
      // write the dataset for the property enum
      status = H5Dwrite(property_enum_dataset,
                        property_enum_id,
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        &property);
      AssertThrow(status >= 0, ExcIO());

      status = H5Dclose(state_enum_dataset);
      AssertThrow(status >= 0, ExcIO());
      status = H5Dclose(property_enum_dataset);
      AssertThrow(status >= 0, ExcIO());
      status = H5Sclose(state_enum_dataspace);
      AssertThrow(status >= 0, ExcIO());
      status = H5Sclose(property_enum_dataspace);
      AssertThrow(status >= 0, ExcIO());
      status = H5Tclose(state_enum_id);
      AssertThrow(status >= 0, ExcIO());
      status = H5Tclose(property_enum_id);
      AssertThrow(status >= 0, ExcIO());
      status = H5Fclose(file_id_reopen);
      AssertThrow(status >= 0, ExcIO());
    }

#  endif
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::load(const std::string &filename)
{
#  ifndef DEAL_II_WITH_HDF5
  (void)filename;
  AssertThrow(false, ExcMessage("HDF5 support is disabled."));
#  else
#    ifdef H5_HAVE_PARALLEL
  // implementation for configurations equipped with a parallel file system
  load_parallel(filename);

#    else
  // implementation for configurations with no parallel file system
  load_serial(filename);
#    endif
#  endif
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::load_serial(const std::string &filename)
{
#  ifndef DEAL_II_WITH_HDF5
  (void)filename;
  Assert(false, ExcInternalError());
#  else

  /*
   * The content of the distributed matrix is copied to a matrix using a 1x1
   * process grid. Therefore, one process has all the data and can write it to a
   * file
   */
  // create a 1xP column grid with P being the number of MPI processes
  const auto one_grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(this->grid->mpi_communicator,
                                                  1,
                                                  1);

  const int                   MB = n_rows, NB = n_columns;
  ScaLAPACKMatrix<NumberType> tmp(n_rows, n_columns, one_grid, MB, NB);

  int state_int    = -1;
  int property_int = -1;

  // the 1x1 grid has only one process and this one reads
  // the content of the matrix from the HDF5 file
  if (tmp.grid->mpi_process_is_active)
    {
      herr_t status;

      // open the file in read-only mode
      hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

      // open the dataset in the file
      hid_t dataset_id = H5Dopen2(file_id, "/matrix", H5P_DEFAULT);

      // check the datatype of the data in the file
      // datatype of source and destination must have the same class
      // see HDF User's Guide: 6.10. Data Transfer: Datatype Conversion and
      // Selection
      hid_t       datatype   = H5Dget_type(dataset_id);
      H5T_class_t t_class_in = H5Tget_class(datatype);
      H5T_class_t t_class    = H5Tget_class(hdf5_type_id(tmp.values.data()));
      AssertThrow(
        t_class_in == t_class,
        ExcMessage(
          "The data type of the matrix to be read does not match the archive"));

      // get dataspace handle
      hid_t dataspace_id = H5Dget_space(dataset_id);
      // get number of dimensions
      const int ndims = H5Sget_simple_extent_ndims(dataspace_id);
      AssertThrow(ndims == 2, ExcIO());
      // get every dimension
      hsize_t dims[2];
      H5Sget_simple_extent_dims(dataspace_id, dims, nullptr);
      AssertThrow(
        static_cast<int>(dims[0]) == n_columns,
        ExcMessage(
          "The number of columns of the matrix does not match the content of the archive"));
      AssertThrow(
        static_cast<int>(dims[1]) == n_rows,
        ExcMessage(
          "The number of rows of the matrix does not match the content of the archive"));

      // read data
      status = H5Dread(dataset_id,
                       hdf5_type_id(tmp.values.data()),
                       H5S_ALL,
                       H5S_ALL,
                       H5P_DEFAULT,
                       tmp.values.data());
      AssertThrow(status >= 0, ExcIO());

      // create HDF5 enum type for LAPACKSupport::State and
      // LAPACKSupport::Property
      hid_t state_enum_id, property_enum_id;
      internal::create_HDF5_state_enum_id(state_enum_id);
      internal::create_HDF5_property_enum_id(property_enum_id);

      // open the datasets for the state and property enum in the file
      hid_t       dataset_state_id = H5Dopen2(file_id, "/state", H5P_DEFAULT);
      hid_t       datatype_state   = H5Dget_type(dataset_state_id);
      H5T_class_t t_class_state    = H5Tget_class(datatype_state);
      AssertThrow(t_class_state == H5T_ENUM, ExcIO());

      hid_t dataset_property_id = H5Dopen2(file_id, "/property", H5P_DEFAULT);
      hid_t datatype_property   = H5Dget_type(dataset_property_id);
      H5T_class_t t_class_property = H5Tget_class(datatype_property);
      AssertThrow(t_class_property == H5T_ENUM, ExcIO());

      // get dataspace handles
      hid_t dataspace_state    = H5Dget_space(dataset_state_id);
      hid_t dataspace_property = H5Dget_space(dataset_property_id);
      // get number of dimensions
      const int ndims_state = H5Sget_simple_extent_ndims(dataspace_state);
      AssertThrow(ndims_state == 1, ExcIO());
      const int ndims_property = H5Sget_simple_extent_ndims(dataspace_property);
      AssertThrow(ndims_property == 1, ExcIO());
      // get every dimension
      hsize_t dims_state[1];
      H5Sget_simple_extent_dims(dataspace_state, dims_state, nullptr);
      AssertThrow(static_cast<int>(dims_state[0]) == 1, ExcIO());
      hsize_t dims_property[1];
      H5Sget_simple_extent_dims(dataspace_property, dims_property, nullptr);
      AssertThrow(static_cast<int>(dims_property[0]) == 1, ExcIO());

      // read data
      status = H5Dread(dataset_state_id,
                       state_enum_id,
                       H5S_ALL,
                       H5S_ALL,
                       H5P_DEFAULT,
                       &tmp.state);
      AssertThrow(status >= 0, ExcIO());
      // To send the state from the root process to the other processes
      // the state enum is casted to an integer, that will be broadcasted and
      // subsequently casted back to the enum type
      state_int = static_cast<int>(tmp.state);

      status = H5Dread(dataset_property_id,
                       property_enum_id,
                       H5S_ALL,
                       H5S_ALL,
                       H5P_DEFAULT,
                       &tmp.property);
      AssertThrow(status >= 0, ExcIO());
      // To send the property from the root process to the other processes
      // the state enum is casted to an integer, that will be broadcasted and
      // subsequently casted back to the enum type
      property_int = static_cast<int>(tmp.property);

      // terminate access to the data spaces
      status = H5Sclose(dataspace_id);
      AssertThrow(status >= 0, ExcIO());
      status = H5Sclose(dataspace_state);
      AssertThrow(status >= 0, ExcIO());
      status = H5Sclose(dataspace_property);
      AssertThrow(status >= 0, ExcIO());

      // release data type handles
      status = H5Tclose(datatype);
      AssertThrow(status >= 0, ExcIO());
      status = H5Tclose(state_enum_id);
      AssertThrow(status >= 0, ExcIO());
      status = H5Tclose(property_enum_id);
      AssertThrow(status >= 0, ExcIO());

      // end access to the data sets and release resources used by them
      status = H5Dclose(dataset_state_id);
      AssertThrow(status >= 0, ExcIO());
      status = H5Dclose(dataset_id);
      AssertThrow(status >= 0, ExcIO());
      status = H5Dclose(dataset_property_id);
      AssertThrow(status >= 0, ExcIO());

      // close the file.
      status = H5Fclose(file_id);
      AssertThrow(status >= 0, ExcIO());
    }
  // so far only the root process has the correct state integer --> broadcasting
  tmp.grid->send_to_inactive(&state_int, 1);
  // so far only the root process has the correct property integer -->
  // broadcasting
  tmp.grid->send_to_inactive(&property_int, 1);

  tmp.state    = static_cast<LAPACKSupport::State>(state_int);
  tmp.property = static_cast<LAPACKSupport::Property>(property_int);

  tmp.copy_to(*this);

#  endif // DEAL_II_WITH_HDF5
}



template <typename NumberType>
void
ScaLAPACKMatrix<NumberType>::load_parallel(const std::string &filename)
{
#  ifndef DEAL_II_WITH_HDF5
  (void)filename;
  Assert(false, ExcInternalError());
#  else
#    ifndef H5_HAVE_PARALLEL
  Assert(false, ExcInternalError());
#    else

  const unsigned int n_mpi_processes(
    Utilities::MPI::n_mpi_processes(this->grid->mpi_communicator));
  MPI_Info info = MPI_INFO_NULL;
  /*
   * The content of the distributed matrix is copied to a matrix using a
   * 1xn_processes process grid. Therefore, the processes hold contiguous chunks
   * of the matrix, which they can write to the file
   */
  // create a 1xP column grid with P being the number of MPI processes
  const auto column_grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(this->grid->mpi_communicator,
                                                  1,
                                                  n_mpi_processes);

  const int MB = n_rows;
  // for the choice of NB see explanation in save_parallel()
  const int NB = std::max(static_cast<int>(std::ceil(
                            static_cast<double>(n_columns) / n_mpi_processes)),
                          column_block_size);

  ScaLAPACKMatrix<NumberType> tmp(n_rows, n_columns, column_grid, MB, NB);

  // get pointer to data held by the process
  NumberType *data = (tmp.values.size() > 0) ? tmp.values.data() : nullptr;

  herr_t status;

  // set up file access property list with parallel I/O access
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  status         = H5Pset_fapl_mpio(plist_id, tmp.grid->mpi_communicator, info);
  AssertThrow(status >= 0, ExcIO());

  // open file collectively in read-only mode and release property list
  // identifier
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist_id);
  status        = H5Pclose(plist_id);
  AssertThrow(status >= 0, ExcIO());

  // open the dataset in the file collectively
  hid_t dataset_id = H5Dopen2(file_id, "/matrix", H5P_DEFAULT);

  // check the datatype of the dataset in the file
  // if the classes of type of the dataset and the matrix do not match abort
  // see HDF User's Guide: 6.10. Data Transfer: Datatype Conversion and
  // Selection
  hid_t       datatype     = hdf5_type_id(data);
  hid_t       datatype_inp = H5Dget_type(dataset_id);
  H5T_class_t t_class_inp  = H5Tget_class(datatype_inp);
  H5T_class_t t_class      = H5Tget_class(datatype);
  AssertThrow(
    t_class_inp == t_class,
    ExcMessage(
      "The data type of the matrix to be read does not match the archive"));

  // get the dimensions of the matrix stored in the file
  // get dataspace handle
  hid_t dataspace_id = H5Dget_space(dataset_id);
  // get number of dimensions
  const int ndims = H5Sget_simple_extent_ndims(dataspace_id);
  AssertThrow(ndims == 2, ExcIO());
  // get every dimension
  hsize_t dims[2];
  status = H5Sget_simple_extent_dims(dataspace_id, dims, nullptr);
  AssertThrow(status >= 0, ExcIO());
  AssertThrow(
    static_cast<int>(dims[0]) == n_columns,
    ExcMessage(
      "The number of columns of the matrix does not match the content of the archive"));
  AssertThrow(
    static_cast<int>(dims[1]) == n_rows,
    ExcMessage(
      "The number of rows of the matrix does not match the content of the archive"));

  // gather the number of local rows and columns from all processes
  std::vector<int> proc_n_local_rows(n_mpi_processes),
    proc_n_local_columns(n_mpi_processes);
  int ierr = MPI_Allgather(&tmp.n_local_rows,
                           1,
                           MPI_INT,
                           proc_n_local_rows.data(),
                           1,
                           MPI_INT,
                           tmp.grid->mpi_communicator);
  AssertThrowMPI(ierr);
  ierr = MPI_Allgather(&tmp.n_local_columns,
                       1,
                       MPI_INT,
                       proc_n_local_columns.data(),
                       1,
                       MPI_INT,
                       tmp.grid->mpi_communicator);
  AssertThrowMPI(ierr);

  const unsigned int my_rank(
    Utilities::MPI::this_mpi_process(tmp.grid->mpi_communicator));

  // hyperslab selection parameters
  // each process defines dataset in memory and writes it to the hyperslab in
  // the file
  hsize_t count[2];
  count[0] = tmp.n_local_columns;
  count[1] = tmp.n_local_rows;

  hsize_t offset[2] = {0};
  for (unsigned int i = 0; i < my_rank; ++i)
    offset[0] += proc_n_local_columns[i];

  // select hyperslab in the file
  status = H5Sselect_hyperslab(
    dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
  AssertThrow(status >= 0, ExcIO());

  // create a memory dataspace independently
  hid_t memspace = H5Screate_simple(2, count, nullptr);

  // read data independently
  status =
    H5Dread(dataset_id, datatype, memspace, dataspace_id, H5P_DEFAULT, data);
  AssertThrow(status >= 0, ExcIO());

  // create HDF5 enum type for LAPACKSupport::State and LAPACKSupport::Property
  hid_t state_enum_id, property_enum_id;
  internal::create_HDF5_state_enum_id(state_enum_id);
  internal::create_HDF5_property_enum_id(property_enum_id);

  // open the datasets for the state and property enum in the file
  hid_t       dataset_state_id = H5Dopen2(file_id, "/state", H5P_DEFAULT);
  hid_t       datatype_state   = H5Dget_type(dataset_state_id);
  H5T_class_t t_class_state    = H5Tget_class(datatype_state);
  AssertThrow(t_class_state == H5T_ENUM, ExcIO());

  hid_t       dataset_property_id = H5Dopen2(file_id, "/property", H5P_DEFAULT);
  hid_t       datatype_property   = H5Dget_type(dataset_property_id);
  H5T_class_t t_class_property    = H5Tget_class(datatype_property);
  AssertThrow(t_class_property == H5T_ENUM, ExcIO());

  // get dataspace handles
  hid_t dataspace_state    = H5Dget_space(dataset_state_id);
  hid_t dataspace_property = H5Dget_space(dataset_property_id);
  // get number of dimensions
  const int ndims_state = H5Sget_simple_extent_ndims(dataspace_state);
  AssertThrow(ndims_state == 1, ExcIO());
  const int ndims_property = H5Sget_simple_extent_ndims(dataspace_property);
  AssertThrow(ndims_property == 1, ExcIO());
  // get every dimension
  hsize_t dims_state[1];
  H5Sget_simple_extent_dims(dataspace_state, dims_state, nullptr);
  AssertThrow(static_cast<int>(dims_state[0]) == 1, ExcIO());
  hsize_t dims_property[1];
  H5Sget_simple_extent_dims(dataspace_property, dims_property, nullptr);
  AssertThrow(static_cast<int>(dims_property[0]) == 1, ExcIO());

  // read data
  status = H5Dread(
    dataset_state_id, state_enum_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmp.state);
  AssertThrow(status >= 0, ExcIO());

  status = H5Dread(dataset_property_id,
                   property_enum_id,
                   H5S_ALL,
                   H5S_ALL,
                   H5P_DEFAULT,
                   &tmp.property);
  AssertThrow(status >= 0, ExcIO());

  // close/release sources
  status = H5Sclose(memspace);
  AssertThrow(status >= 0, ExcIO());
  status = H5Dclose(dataset_id);
  AssertThrow(status >= 0, ExcIO());
  status = H5Dclose(dataset_state_id);
  AssertThrow(status >= 0, ExcIO());
  status = H5Dclose(dataset_property_id);
  AssertThrow(status >= 0, ExcIO());
  status = H5Sclose(dataspace_id);
  AssertThrow(status >= 0, ExcIO());
  status = H5Sclose(dataspace_state);
  AssertThrow(status >= 0, ExcIO());
  status = H5Sclose(dataspace_property);
  AssertThrow(status >= 0, ExcIO());
  // status = H5Tclose(datatype);
  // AssertThrow(status >= 0, ExcIO());
  status = H5Tclose(state_enum_id);
  AssertThrow(status >= 0, ExcIO());
  status = H5Tclose(property_enum_id);
  AssertThrow(status >= 0, ExcIO());
  status = H5Fclose(file_id);
  AssertThrow(status >= 0, ExcIO());

  // copying the distributed matrices
  tmp.copy_to(*this);

#    endif // H5_HAVE_PARALLEL
#  endif   // DEAL_II_WITH_HDF5
}



namespace internal
{
  namespace
  {
    template <typename NumberType>
    void
    scale_columns(ScaLAPACKMatrix<NumberType> &      matrix,
                  const ArrayView<const NumberType> &factors)
    {
      Assert(matrix.n() == factors.size(),
             ExcDimensionMismatch(matrix.n(), factors.size()));

      for (unsigned int i = 0; i < matrix.local_n(); ++i)
        {
          const NumberType s = factors[matrix.global_column(i)];

          for (unsigned int j = 0; j < matrix.local_m(); ++j)
            matrix.local_el(j, i) *= s;
        }
    }

    template <typename NumberType>
    void
    scale_rows(ScaLAPACKMatrix<NumberType> &      matrix,
               const ArrayView<const NumberType> &factors)
    {
      Assert(matrix.m() == factors.size(),
             ExcDimensionMismatch(matrix.m(), factors.size()));

      for (unsigned int i = 0; i < matrix.local_m(); ++i)
        {
          const NumberType s = factors[matrix.global_row(i)];

          for (unsigned int j = 0; j < matrix.local_n(); ++j)
            matrix.local_el(i, j) *= s;
        }
    }

  } // namespace
} // namespace internal



template <typename NumberType>
template <class InputVector>
void
ScaLAPACKMatrix<NumberType>::scale_columns(const InputVector &factors)
{
  if (this->grid->mpi_process_is_active)
    internal::scale_columns(*this, make_array_view(factors));
}



template <typename NumberType>
template <class InputVector>
void
ScaLAPACKMatrix<NumberType>::scale_rows(const InputVector &factors)
{
  if (this->grid->mpi_process_is_active)
    internal::scale_rows(*this, make_array_view(factors));
}



// instantiations
#  include "scalapack.inst"


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SCALAPACK
