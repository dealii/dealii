// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/sparsity_pattern.h>
#  include <deal.II/lac/compressed_sparsity_pattern.h>
#  include <deal.II/lac/compressed_simple_sparsity_pattern.h>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  namespace MPI
  {

    SparseMatrix::SparseMatrix ()
    {
      // just like for vectors: since we
      // create an empty matrix, we can as
      // well make it sequential
      const int m=0, n=0, n_nonzero_per_row=0;
      const int ierr
        = MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, n_nonzero_per_row,
                          0, &matrix);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
    }


    SparseMatrix::~SparseMatrix ()
    {
      int ierr;

#if DEAL_II_PETSC_VERSION_LT(3,2,0)
      ierr = MatDestroy (matrix);
#else
      ierr = MatDestroy (&matrix);
#endif
      AssertThrow (ierr == 0, ExcPETScError(ierr));
    }

    SparseMatrix::SparseMatrix (const MPI_Comm  &communicator,
                                const size_type  m,
                                const size_type  n,
                                const size_type  local_rows,
                                const size_type  local_columns,
                                const size_type  n_nonzero_per_row,
                                const bool       is_symmetric)
      :
      communicator (communicator)
    {
      do_reinit (m, n, local_rows, local_columns,
                 n_nonzero_per_row, is_symmetric);
    }



    SparseMatrix::SparseMatrix (const MPI_Comm               &communicator,
                                const size_type               m,
                                const size_type               n,
                                const size_type               local_rows,
                                const size_type               local_columns,
                                const std::vector<size_type> &row_lengths,
                                const bool                    is_symmetric)
      :
      communicator (communicator)
    {
      do_reinit (m, n, local_rows, local_columns,
                 row_lengths, is_symmetric);
    }



    template <typename SparsityType>
    SparseMatrix::
    SparseMatrix (const MPI_Comm               &communicator,
                  const SparsityType           &sparsity_pattern,
                  const std::vector<size_type> &local_rows_per_process,
                  const std::vector<size_type> &local_columns_per_process,
                  const unsigned int            this_process,
                  const bool                    preset_nonzero_locations)
      :
      communicator (communicator)
    {
      do_reinit (sparsity_pattern, local_rows_per_process,
                 local_columns_per_process, this_process,
                 preset_nonzero_locations);
    }



    SparseMatrix &
    SparseMatrix::operator = (const value_type d)
    {
      MatrixBase::operator = (d);
      return *this;
    }

    void
    SparseMatrix::copy_from (const SparseMatrix &other)
    {
      if (&other == this)
        return;

      this->communicator = other.communicator;

      int ierr = MatCopy(other.matrix, matrix, SAME_NONZERO_PATTERN);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
    }

    void
    SparseMatrix::reinit (const MPI_Comm  &communicator,
                          const size_type  m,
                          const size_type  n,
                          const size_type  local_rows,
                          const size_type  local_columns,
                          const size_type  n_nonzero_per_row,
                          const bool       is_symmetric)
    {
      this->communicator = communicator;

      // get rid of old matrix and generate a
      // new one
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
      const int ierr = MatDestroy (matrix);
#else
      const int ierr = MatDestroy (&matrix);
#endif
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      do_reinit (m, n, local_rows, local_columns,
                 n_nonzero_per_row, is_symmetric);
    }



    void
    SparseMatrix::reinit (const MPI_Comm               &communicator,
                          const size_type               m,
                          const size_type               n,
                          const size_type               local_rows,
                          const size_type               local_columns,
                          const std::vector<size_type> &row_lengths,
                          const bool                    is_symmetric)
    {
      this->communicator = communicator;

      // get rid of old matrix and generate a
      // new one
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
      const int ierr = MatDestroy (matrix);
#else
      const int ierr = MatDestroy (&matrix);
#endif
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      do_reinit (m, n, local_rows, local_columns, row_lengths, is_symmetric);
    }



    template <typename SparsityType>
    void
    SparseMatrix::
    reinit (const MPI_Comm               &communicator,
            const SparsityType           &sparsity_pattern,
            const std::vector<size_type> &local_rows_per_process,
            const std::vector<size_type> &local_columns_per_process,
            const unsigned int            this_process,
            const bool                    preset_nonzero_locations)
    {
      this->communicator = communicator;

      // get rid of old matrix and generate a
      // new one
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
      const int ierr = MatDestroy (matrix);
#else
      const int ierr = MatDestroy (&matrix);
#endif
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      do_reinit (sparsity_pattern, local_rows_per_process,
                 local_columns_per_process, this_process,
                 preset_nonzero_locations);
    }

    template <typename SparsityType>
    void
    SparseMatrix::
    reinit (const IndexSet &local_rows,
            const IndexSet &local_columns,
            const SparsityType &sparsity_pattern,
            const MPI_Comm &communicator)
    {
      this->communicator = communicator;

      // get rid of old matrix and generate a
      // new one
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
      const int ierr = MatDestroy (matrix);
#else
      const int ierr = MatDestroy (&matrix);
#endif
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      do_reinit (local_rows, local_columns, sparsity_pattern);
    }

    void
    SparseMatrix::do_reinit (const size_type m,
                             const size_type n,
                             const size_type local_rows,
                             const size_type local_columns,
                             const size_type n_nonzero_per_row,
                             const bool      is_symmetric)
    {
      Assert (local_rows <= m, ExcLocalRowsTooLarge (local_rows, m));

      // use the call sequence indicating only
      // a maximal number of elements per row
      // for all rows globally
      int ierr;

#if DEAL_II_PETSC_VERSION_LT(3,3,0)
      ierr
        = MatCreateMPIAIJ (communicator,
                           local_rows, local_columns,
                           m, n,
                           n_nonzero_per_row, 0, 0, 0,
                           &matrix);
#else
      ierr
        = MatCreateAIJ (communicator,
                        local_rows, local_columns,
                        m, n,
                        n_nonzero_per_row, 0, 0, 0,
                        &matrix);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      ierr = MatSetOption (matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
#endif
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      // set symmetric flag, if so requested
      if (is_symmetric == true)
        {
#if DEAL_II_PETSC_VERSION_LT(3,0,0)
          const int ierr
            = MatSetOption (matrix, MAT_SYMMETRIC);
#else
          const int ierr
            = MatSetOption (matrix, MAT_SYMMETRIC, PETSC_TRUE);
#endif

          AssertThrow (ierr == 0, ExcPETScError(ierr));
        }
    }



    void
    SparseMatrix::do_reinit (const size_type m,
                             const size_type n,
                             const size_type local_rows,
                             const size_type local_columns,
                             const std::vector<size_type> &row_lengths,
                             const bool      is_symmetric)
    {
      Assert (local_rows <= m, ExcLocalRowsTooLarge (local_rows, m));

      Assert (row_lengths.size() == m,
              ExcDimensionMismatch (row_lengths.size(), m));

      // For the case that
      // local_columns is smaller
      // than one of the row lengths
      // MatCreateMPIAIJ throws an
      // error. In this case use a
      // PETScWrappers::SparseMatrix
      for (size_type i=0; i<row_lengths.size(); ++i)
        Assert(row_lengths[i]<=local_columns,
               ExcIndexRange(row_lengths[i], 1, local_columns+1));

      // use the call sequence indicating a
      // maximal number of elements for each
      // row individually. annoyingly, we
      // always use unsigned ints for cases
      // like this, while PETSc wants to see
      // signed integers. so we have to
      // convert, unless we want to play dirty
      // tricks with conversions of pointers
      const std::vector<PetscInt> int_row_lengths (row_lengths.begin(),
                                                   row_lengths.end());

//TODO: There must be a significantly better way to provide information about the off-diagonal blocks of the matrix. this way, petsc keeps allocating tiny chunks of memory, and gets completely hung up over this

      int ierr;

#if DEAL_II_PETSC_VERSION_LT(3,3,0)
      ierr
        = MatCreateMPIAIJ (communicator,
                           local_rows, local_columns,
                           m, n,
                           0, &int_row_lengths[0], 0, 0,
                           &matrix);
#else
      ierr
        = MatCreateAIJ (communicator,
                        local_rows, local_columns,
                        m, n,
                        0, &int_row_lengths[0], 0, 0,
                        &matrix);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

//TODO: Sometimes the actual number of nonzero entries allocated is greater than the number of nonzero entries, which petsc will complain about unless explicitly disabled with MatSetOption. There is probably a way to prevent a different number nonzero elements being allocated in the first place. (See also previous TODO).

      ierr = MatSetOption (matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
#endif
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      // set symmetric flag, if so requested
      if (is_symmetric == true)
        {
#if DEAL_II_PETSC_VERSION_LT(3,0,0)
          const int ierr
            = MatSetOption (matrix, MAT_SYMMETRIC);
#else
          const int ierr
            = MatSetOption (matrix, MAT_SYMMETRIC, PETSC_TRUE);
#endif

          AssertThrow (ierr == 0, ExcPETScError(ierr));
        }
    }


    template <typename SparsityType>
    void
    SparseMatrix::
    do_reinit (const IndexSet &local_rows,
               const IndexSet &local_columns,
               const SparsityType         &sparsity_pattern)
    {
      Assert(sparsity_pattern.n_rows()==local_rows.size(),
             ExcMessage("SparsityPattern and IndexSet have different number of rows"));
      Assert(sparsity_pattern.n_cols()==local_columns.size(),
             ExcMessage("SparsityPattern and IndexSet have different number of columns"));
      Assert(local_rows.is_contiguous() && local_columns.is_contiguous(),
             ExcMessage("PETSc only supports contiguous row/column ranges"));




      // create the matrix. We do not set row length but set the
      // correct SparsityPattern later.
      int ierr;

      ierr = MatCreate(communicator,&matrix);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      ierr = MatSetSizes(matrix,
                         local_rows.n_elements(),
                         local_columns.n_elements(),
                         sparsity_pattern.n_rows(),
                         sparsity_pattern.n_cols());
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      ierr = MatSetType(matrix,MATMPIAIJ);
      AssertThrow (ierr == 0, ExcPETScError(ierr));


      // next preset the exact given matrix
      // entries with zeros. this doesn't avoid any
      // memory allocations, but it at least
      // avoids some searches later on. the
      // key here is that we can use the
      // matrix set routines that set an
      // entire row at once, not a single
      // entry at a time
      //
      // for the usefulness of this option
      // read the documentation of this
      // class.
      //if (preset_nonzero_locations == true)
      if (local_rows.n_elements()>0)
        {
          Assert(local_columns.n_elements()>0, ExcInternalError());
          // MatMPIAIJSetPreallocationCSR
          // can be used to allocate the sparsity
          // pattern of a matrix

          const PetscInt local_row_start = local_rows.nth_index_in_set(0);
          const PetscInt
          local_row_end = local_row_start + local_rows.n_elements();


          // first set up the column number
          // array for the rows to be stored
          // on the local processor. have one
          // dummy entry at the end to make
          // sure petsc doesn't read past the
          // end
          std::vector<PetscInt>

          rowstart_in_window (local_row_end - local_row_start + 1, 0),
                             colnums_in_window;
          {
            unsigned int n_cols = 0;
            for (PetscInt i=local_row_start; i<local_row_end; ++i)
              {
                const PetscInt row_length = sparsity_pattern.row_length(i);
                rowstart_in_window[i+1-local_row_start]
                  = rowstart_in_window[i-local_row_start] + row_length;
                n_cols += row_length;
              }
            colnums_in_window.resize (n_cols+1, -1);
          }

          // now copy over the information
          // from the sparsity pattern.
          {
            PetscInt *ptr = & colnums_in_window[0];

            for (PetscInt i=local_row_start; i<local_row_end; ++i)
              {
                typename SparsityType::row_iterator
                row_start = sparsity_pattern.row_begin(i),
                row_end = sparsity_pattern.row_end(i);

                std::copy(row_start, row_end, ptr);
                ptr += row_end - row_start;
              }
          }


          // then call the petsc function
          // that summarily allocates these
          // entries:
          MatMPIAIJSetPreallocationCSR (matrix,
                                        &rowstart_in_window[0],
                                        &colnums_in_window[0],
                                        0);
        }
      else
        {
          PetscInt i=0;
          MatMPIAIJSetPreallocationCSR (matrix,
                                        &i,
                                        &i,
                                        0);


        }
      compress (dealii::VectorOperation::insert);

      {

        // Tell PETSc that we are not
        // planning on adding new entries
        // to the matrix. Generate errors
        // in debug mode.
        int ierr;
#if DEAL_II_PETSC_VERSION_LT(3,0,0)
#ifdef DEBUG
        ierr = MatSetOption (matrix, MAT_NEW_NONZERO_LOCATION_ERR);
        AssertThrow (ierr == 0, ExcPETScError(ierr));
#else
        ierr = MatSetOption (matrix, MAT_NO_NEW_NONZERO_LOCATIONS);
        AssertThrow (ierr == 0, ExcPETScError(ierr));
#endif
#else
#ifdef DEBUG
        ierr = MatSetOption (matrix, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
        AssertThrow (ierr == 0, ExcPETScError(ierr));
#else
        ierr = MatSetOption (matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
        AssertThrow (ierr == 0, ExcPETScError(ierr));
#endif
#endif

        // Tell PETSc to keep the
        // SparsityPattern entries even if
        // we delete a row with
        // clear_rows() which calls
        // MatZeroRows(). Otherwise one can
        // not write into that row
        // afterwards.
#if DEAL_II_PETSC_VERSION_LT(3,0,0)
        ierr = MatSetOption (matrix, MAT_KEEP_ZEROED_ROWS);
        AssertThrow (ierr == 0, ExcPETScError(ierr));
#elif DEAL_II_PETSC_VERSION_LT(3,1,0)
        ierr = MatSetOption (matrix, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE);
        AssertThrow (ierr == 0, ExcPETScError(ierr));
#else
        ierr = MatSetOption (matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
        AssertThrow (ierr == 0, ExcPETScError(ierr));
#endif

      }

    }


    template <typename SparsityType>
    void
    SparseMatrix::
    do_reinit (const SparsityType           &sparsity_pattern,
               const std::vector<size_type> &local_rows_per_process,
               const std::vector<size_type> &local_columns_per_process,
               const unsigned int            this_process,
               const bool                    preset_nonzero_locations)
    {
      Assert (local_rows_per_process.size() == local_columns_per_process.size(),
              ExcDimensionMismatch (local_rows_per_process.size(),
                                    local_columns_per_process.size()));
      Assert (this_process < local_rows_per_process.size(),
              ExcInternalError());
      // for each row that we own locally, we
      // have to count how many of the
      // entries in the sparsity pattern lie
      // in the column area we have locally,
      // and how many arent. for this, we
      // first have to know which areas are
      // ours
      size_type local_row_start = 0;
      size_type local_col_start = 0;
      for (unsigned int p=0; p<this_process; ++p)
        {
          local_row_start += local_rows_per_process[p];
          local_col_start += local_columns_per_process[p];
        }
      const size_type
      local_row_end = local_row_start + local_rows_per_process[this_process];

#if DEAL_II_PETSC_VERSION_LT(2,3,3)
      //old version to create the matrix, we
      //can skip calculating the row length
      //at least starting from 2.3.3 (tested,
      //see below)

      const size_type
      local_col_end = local_col_start + local_columns_per_process[this_process];

      // then count the elements in- and
      // out-of-window for the rows we own
      std::vector<PetscInt>

      row_lengths_in_window (local_row_end - local_row_start),
                            row_lengths_out_of_window (local_row_end - local_row_start);
      for (size_type row = local_row_start; row<local_row_end; ++row)
        for (size_type c=0; c<sparsity_pattern.row_length(row); ++c)
          {
            const size_type column = sparsity_pattern.column_number(row,c);

            if ((column >= local_col_start) &&
                (column < local_col_end))
              ++row_lengths_in_window[row-local_row_start];
            else
              ++row_lengths_out_of_window[row-local_row_start];
          }


      // create the matrix. completely
      // confusingly, PETSc wants us to pass
      // arrays for the local number of
      // elements that starts with zero for
      // the first _local_ row, i.e. it
      // doesn't index into an array for
      // _all_ rows.
      const int ierr
        = MatCreateMPIAIJ(communicator,
                          local_rows_per_process[this_process],
                          local_columns_per_process[this_process],
                          sparsity_pattern.n_rows(),
                          sparsity_pattern.n_cols(),
                          0, &row_lengths_in_window[0],
                          0, &row_lengths_out_of_window[0],
                          &matrix);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

#else //PETSC_VERSION>=2.3.3
      // create the matrix. We
      // do not set row length but set the
      // correct SparsityPattern later.
      int ierr;

      ierr = MatCreate(communicator,&matrix);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      ierr = MatSetSizes(matrix,
                         local_rows_per_process[this_process],
                         local_columns_per_process[this_process],
                         sparsity_pattern.n_rows(),
                         sparsity_pattern.n_cols());
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      ierr = MatSetType(matrix,MATMPIAIJ);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
#endif


      // next preset the exact given matrix
      // entries with zeros, if the user
      // requested so. this doesn't avoid any
      // memory allocations, but it at least
      // avoids some searches later on. the
      // key here is that we can use the
      // matrix set routines that set an
      // entire row at once, not a single
      // entry at a time
      //
      // for the usefulness of this option
      // read the documentation of this
      // class.
      if (preset_nonzero_locations == true)
        {
          // MatMPIAIJSetPreallocationCSR
          // can be used to allocate the sparsity
          // pattern of a matrix if it is already
          // available:

          // first set up the column number
          // array for the rows to be stored
          // on the local processor. have one
          // dummy entry at the end to make
          // sure petsc doesn't read past the
          // end
          std::vector<PetscInt>

          rowstart_in_window (local_row_end - local_row_start + 1, 0),
                             colnums_in_window;
          {
            size_type n_cols = 0;
            for (size_type i=local_row_start; i<local_row_end; ++i)
              {
                const size_type row_length = sparsity_pattern.row_length(i);
                rowstart_in_window[i+1-local_row_start]
                  = rowstart_in_window[i-local_row_start] + row_length;
                n_cols += row_length;
              }
            colnums_in_window.resize (n_cols+1, -1);
          }

          // now copy over the information
          // from the sparsity pattern.
          {
            PetscInt *ptr = & colnums_in_window[0];

            for (size_type i=local_row_start; i<local_row_end; ++i)
              {
                typename SparsityType::row_iterator
                row_start = sparsity_pattern.row_begin(i),
                row_end = sparsity_pattern.row_end(i);

                std::copy(row_start, row_end, ptr);
                ptr += row_end - row_start;
              }
          }


          // then call the petsc function
          // that summarily allocates these
          // entries:
          MatMPIAIJSetPreallocationCSR (matrix,
                                        &rowstart_in_window[0],
                                        &colnums_in_window[0],
                                        0);

#if DEAL_II_PETSC_VERSION_LT(2,3,3)
          // this is only needed for old
          // PETSc versions:

          // for some reason, it does not
          // seem to be possible to force
          // actual allocation of actual
          // entries by using the last
          // arguments to the call above. if
          // we don't initialize the entries
          // like in the following loop, then
          // the program is unbearably slow
          // because elements are allocated
          // and accessed in random order,
          // which is not what PETSc likes
          //
          // note that we actually have to
          // set the entries to something
          // non-zero! do the allocation one
          // row at a time
          {
            const std::vector<PetscScalar>
            values (sparsity_pattern.max_entries_per_row(),
                    1.);

            for (size_type i=local_row_start; i<local_row_end; ++i)
              {
                PetscInt petsc_i = i;
                MatSetValues (matrix, 1, &petsc_i,
                              sparsity_pattern.row_length(i),
                              &colnums_in_window[rowstart_in_window[i-local_row_start]],
                              &values[0], INSERT_VALUES);
              }
          }

          compress ();

          // set the dummy entries set above
          // back to zero
          *this = 0;
#endif // version <=2.3.3
          compress ();


          // Tell PETSc that we are not
          // planning on adding new entries
          // to the matrix. Generate errors
          // in debug mode.
          int ierr;
#if DEAL_II_PETSC_VERSION_LT(3,0,0)
#ifdef DEBUG
          ierr = MatSetOption (matrix, MAT_NEW_NONZERO_LOCATION_ERR);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
#else
          ierr = MatSetOption (matrix, MAT_NO_NEW_NONZERO_LOCATIONS);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
#endif
#else
#ifdef DEBUG
          ierr = MatSetOption (matrix, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
#else
          ierr = MatSetOption (matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
#endif
#endif

          // Tell PETSc to keep the
          // SparsityPattern entries even if
          // we delete a row with
          // clear_rows() which calls
          // MatZeroRows(). Otherwise one can
          // not write into that row
          // afterwards.
#if DEAL_II_PETSC_VERSION_LT(3,0,0)
          ierr = MatSetOption (matrix, MAT_KEEP_ZEROED_ROWS);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
#elif DEAL_II_PETSC_VERSION_LT(3,1,0)
          ierr = MatSetOption (matrix, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
#else
          ierr = MatSetOption (matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
#endif

        }
    }

    // explicit instantiations
    //
    template
    SparseMatrix::SparseMatrix (const MPI_Comm &,
                                const SparsityPattern &,
                                const std::vector<size_type> &,
                                const std::vector<size_type> &,
                                const unsigned int,
                                const bool);
    template
    SparseMatrix::SparseMatrix (const MPI_Comm &,
                                const CompressedSparsityPattern &,
                                const std::vector<size_type> &,
                                const std::vector<size_type> &,
                                const unsigned int,
                                const bool);
    template
    SparseMatrix::SparseMatrix (const MPI_Comm &,
                                const CompressedSimpleSparsityPattern &,
                                const std::vector<size_type> &,
                                const std::vector<size_type> &,
                                const unsigned int,
                                const bool);

    template void
    SparseMatrix::reinit (const MPI_Comm &,
                          const SparsityPattern &,
                          const std::vector<size_type> &,
                          const std::vector<size_type> &,
                          const unsigned int,
                          const bool);
    template void
    SparseMatrix::reinit (const MPI_Comm &,
                          const CompressedSparsityPattern &,
                          const std::vector<size_type> &,
                          const std::vector<size_type> &,
                          const unsigned int,
                          const bool);
    template void
    SparseMatrix::reinit (const MPI_Comm &,
                          const CompressedSimpleSparsityPattern &,
                          const std::vector<size_type> &,
                          const std::vector<size_type> &,
                          const unsigned int,
                          const bool);

    template void
    SparseMatrix::
    reinit (const IndexSet &,
            const IndexSet &,
            const CompressedSimpleSparsityPattern &,
            const MPI_Comm &);

    template void
    SparseMatrix::do_reinit (const SparsityPattern &,
                             const std::vector<size_type> &,
                             const std::vector<size_type> &,
                             const unsigned int ,
                             const bool);
    template void
    SparseMatrix::do_reinit (const CompressedSparsityPattern &,
                             const std::vector<size_type> &,
                             const std::vector<size_type> &,
                             const unsigned int ,
                             const bool);
    template void
    SparseMatrix::do_reinit (const CompressedSimpleSparsityPattern &,
                             const std::vector<size_type> &,
                             const std::vector<size_type> &,
                             const unsigned int ,
                             const bool);

    template void
    SparseMatrix::
    do_reinit (const IndexSet &,
               const IndexSet &,
               const CompressedSimpleSparsityPattern &);


    PetscScalar
    SparseMatrix::matrix_norm_square (const Vector &v) const
    {
      Vector tmp (v);
      vmult (tmp, v);
      return tmp*v;
    }

    PetscScalar
    SparseMatrix::matrix_scalar_product (const Vector &u,
                                         const Vector &v) const
    {
      Vector tmp (v);
      vmult (tmp, v);
      return u*tmp;
    }

  }
}


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
