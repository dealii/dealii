// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2016 by the deal.II authors
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

#  include <deal.II/base/mpi.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_compatibility.h>
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/sparsity_pattern.h>
#  include <deal.II/lac/dynamic_sparsity_pattern.h>

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
      const PetscErrorCode ierr
        = MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, n_nonzero_per_row,
                          0, &matrix);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
    }


    SparseMatrix::~SparseMatrix ()
    {
      destroy_matrix (matrix);
    }

    SparseMatrix::SparseMatrix (const MPI_Comm  &communicator,
                                const size_type  m,
                                const size_type  n,
                                const size_type  local_rows,
                                const size_type  local_columns,
                                const size_type  n_nonzero_per_row,
                                const bool       is_symmetric,
                                const size_type  n_offdiag_nonzero_per_row)
      :
      communicator (communicator)
    {
      do_reinit (m, n, local_rows, local_columns,
                 n_nonzero_per_row, is_symmetric,
                 n_offdiag_nonzero_per_row);
    }



    SparseMatrix::SparseMatrix (const MPI_Comm               &communicator,
                                const size_type               m,
                                const size_type               n,
                                const size_type               local_rows,
                                const size_type               local_columns,
                                const std::vector<size_type> &row_lengths,
                                const bool                    is_symmetric,
                                const std::vector<size_type> &offdiag_row_lengths)
      :
      communicator (communicator)
    {
      do_reinit (m, n, local_rows, local_columns,
                 row_lengths, is_symmetric, offdiag_row_lengths);
    }



    template <typename SparsityPatternType>
    SparseMatrix::
    SparseMatrix (const MPI_Comm               &communicator,
                  const SparsityPatternType    &sparsity_pattern,
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


    void
    SparseMatrix::
    reinit (const SparseMatrix &other)
    {
      if (&other == this)
        return;

      this->communicator = other.communicator;

      PetscErrorCode ierr = destroy_matrix (matrix);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      ierr = MatDuplicate (other.matrix, MAT_DO_NOT_COPY_VALUES, &matrix);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
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

      const PetscErrorCode ierr = MatCopy(other.matrix, matrix,
                                          SAME_NONZERO_PATTERN);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
    }

    void
    SparseMatrix::reinit (const MPI_Comm  &communicator,
                          const size_type  m,
                          const size_type  n,
                          const size_type  local_rows,
                          const size_type  local_columns,
                          const size_type  n_nonzero_per_row,
                          const bool       is_symmetric,
                          const size_type  n_offdiag_nonzero_per_row)
    {
      this->communicator = communicator;

      // get rid of old matrix and generate a new one
      const PetscErrorCode ierr = destroy_matrix (matrix);
      AssertThrow (ierr == 0, ExcPETScError (ierr));

      do_reinit (m, n, local_rows, local_columns,
                 n_nonzero_per_row, is_symmetric,
                 n_offdiag_nonzero_per_row);
    }



    void
    SparseMatrix::reinit (const MPI_Comm               &communicator,
                          const size_type               m,
                          const size_type               n,
                          const size_type               local_rows,
                          const size_type               local_columns,
                          const std::vector<size_type> &row_lengths,
                          const bool                    is_symmetric,
                          const std::vector<size_type> &offdiag_row_lengths)
    {
      this->communicator = communicator;

      // get rid of old matrix and generate a
      // new one
      const PetscErrorCode ierr = destroy_matrix (matrix);
      AssertThrow (ierr == 0, ExcPETScError (ierr));

      do_reinit (m, n, local_rows, local_columns,
                 row_lengths, is_symmetric, offdiag_row_lengths);
    }



    template <typename SparsityPatternType>
    void
    SparseMatrix::
    reinit (const MPI_Comm               &communicator,
            const SparsityPatternType    &sparsity_pattern,
            const std::vector<size_type> &local_rows_per_process,
            const std::vector<size_type> &local_columns_per_process,
            const unsigned int            this_process,
            const bool                    preset_nonzero_locations)
    {
      this->communicator = communicator;

      // get rid of old matrix and generate a new one
      destroy_matrix (matrix);
      const PetscErrorCode ierr = destroy_matrix (matrix);
      AssertThrow (ierr == 0, ExcPETScError (ierr));


      do_reinit (sparsity_pattern, local_rows_per_process,
                 local_columns_per_process, this_process,
                 preset_nonzero_locations);
    }

    template <typename SparsityPatternType>
    void
    SparseMatrix::
    reinit (const IndexSet            &local_rows,
            const IndexSet            &local_columns,
            const SparsityPatternType &sparsity_pattern,
            const MPI_Comm            &communicator)
    {
      this->communicator = communicator;

      // get rid of old matrix and generate a new one
      const PetscErrorCode ierr = destroy_matrix (matrix);
      AssertThrow(ierr == 0, ExcPETScError (ierr));

      do_reinit (local_rows, local_columns, sparsity_pattern);
    }

    void
    SparseMatrix::do_reinit (const size_type m,
                             const size_type n,
                             const size_type local_rows,
                             const size_type local_columns,
                             const size_type n_nonzero_per_row,
                             const bool      is_symmetric,
                             const size_type n_offdiag_nonzero_per_row)
    {
      Assert (local_rows <= m, ExcLocalRowsTooLarge (local_rows, m));

      // use the call sequence indicating only
      // a maximal number of elements per row
      // for all rows globally
#if DEAL_II_PETSC_VERSION_LT(3,3,0)
      const PetscErrorCode ierr = MatCreateMPIAIJ
                                  (communicator,
                                   local_rows, local_columns,
                                   m, n,
                                   n_nonzero_per_row, 0,
                                   n_offdiag_nonzero_per_row, 0,
                                   &matrix);
#else
      const PetscErrorCode ierr = MatCreateAIJ
                                  (communicator,
                                   local_rows, local_columns,
                                   m, n,
                                   n_nonzero_per_row, 0,
                                   n_offdiag_nonzero_per_row, 0,
                                   &matrix);
      set_matrix_option (matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
#endif
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      // set symmetric flag, if so requested
      if (is_symmetric == true)
        {
          set_matrix_option (matrix, MAT_SYMMETRIC, PETSC_TRUE);
        }
    }



    void
    SparseMatrix::do_reinit (const size_type m,
                             const size_type n,
                             const size_type local_rows,
                             const size_type local_columns,
                             const std::vector<size_type> &row_lengths,
                             const bool      is_symmetric,
                             const std::vector<size_type> &offdiag_row_lengths)
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
      const std::vector<PetscInt> int_offdiag_row_lengths (offdiag_row_lengths.begin(),
                                                           offdiag_row_lengths.end());

//TODO: There must be a significantly better way to provide information about the off-diagonal blocks of the matrix. this way, petsc keeps allocating tiny chunks of memory, and gets completely hung up over this
#if DEAL_II_PETSC_VERSION_LT(3,3,0)
      const PetscErrorCode ierr = MatCreateMPIAIJ
                                  (communicator,
                                   local_rows, local_columns,
                                   m, n,
                                   0, &int_row_lengths[0],
                                   0, offdiag_row_lengths.size() ? &int_offdiag_row_lengths[0] : 0,
                                   &matrix);
#else
      const PetscErrorCode ierr = MatCreateAIJ
                                  (communicator,
                                   local_rows, local_columns,
                                   m, n,
                                   0, &int_row_lengths[0],
                                   0,
                                   offdiag_row_lengths.size() ? &int_offdiag_row_lengths[0] : 0,
                                   &matrix);

//TODO: Sometimes the actual number of nonzero entries allocated is greater than the number of nonzero entries, which petsc will complain about unless explicitly disabled with MatSetOption. There is probably a way to prevent a different number nonzero elements being allocated in the first place. (See also previous TODO).
      set_matrix_option (matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
#endif
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      // set symmetric flag, if so requested
      if (is_symmetric == true)
        {
          set_matrix_option (matrix, MAT_SYMMETRIC, PETSC_TRUE);
        }
    }


    template <typename SparsityPatternType>
    void
    SparseMatrix::
    do_reinit (const IndexSet            &local_rows,
               const IndexSet            &local_columns,
               const SparsityPatternType &sparsity_pattern)
    {
      Assert(sparsity_pattern.n_rows()==local_rows.size(),
             ExcMessage("SparsityPattern and IndexSet have different number of rows"));
      Assert(sparsity_pattern.n_cols()==local_columns.size(),
             ExcMessage("SparsityPattern and IndexSet have different number of columns"));
      Assert(local_rows.is_contiguous() && local_columns.is_contiguous(),
             ExcMessage("PETSc only supports contiguous row/column ranges"));
      Assert(local_rows.is_ascending_and_one_to_one(communicator), ExcNotImplemented());

#ifdef DEBUG
      {
        // check indexsets
        types::global_dof_index row_owners = Utilities::MPI::sum(local_rows.n_elements(), communicator);
        types::global_dof_index col_owners = Utilities::MPI::sum(local_columns.n_elements(), communicator);
        Assert(row_owners == sparsity_pattern.n_rows(),
               ExcMessage(std::string("Each row has to be owned by exactly one owner (n_rows()=")
                          + Utilities::to_string(sparsity_pattern.n_rows())
                          + " but sum(local_rows.n_elements())="
                          + Utilities::to_string(row_owners)
                          + ")"));
        Assert(col_owners == sparsity_pattern.n_cols(),
               ExcMessage(std::string("Each column has to be owned by exactly one owner (n_cols()=")
                          + Utilities::to_string(sparsity_pattern.n_cols())
                          + " but sum(local_columns.n_elements())="
                          + Utilities::to_string(col_owners)
                          + ")"));
      }
#endif


      // create the matrix. We do not set row length but set the
      // correct SparsityPattern later.
      PetscErrorCode ierr = MatCreate(communicator,&matrix);
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
              for (typename SparsityPatternType::iterator p=sparsity_pattern.begin(i);
                   p != sparsity_pattern.end(i); ++p, ++ptr)
                *ptr = p->column();
          }


          // then call the petsc function
          // that summarily allocates these
          // entries:
          ierr = MatMPIAIJSetPreallocationCSR (matrix,
                                               &rowstart_in_window[0],
                                               &colnums_in_window[0],
                                               0);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
        }
      else
        {
          PetscInt i=0;
          ierr = MatMPIAIJSetPreallocationCSR (matrix,
                                               &i,
                                               &i,
                                               0);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
        }
      compress (dealii::VectorOperation::insert);

      {
        close_matrix (matrix);
        set_keep_zero_rows (matrix);
      }
    }


    template <typename SparsityPatternType>
    void
    SparseMatrix::
    do_reinit (const SparsityPatternType    &sparsity_pattern,
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
      assert_is_compressed ();

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
      const PetscErrorCode ierr = MatCreateMPIAIJ(communicator,
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
      PetscErrorCode ierr = MatCreate(communicator,&matrix);
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
              for (typename SparsityPatternType::iterator p=sparsity_pattern.begin(i);
                   p != sparsity_pattern.end(i); ++p, ++ptr)
                *ptr = p->column();
          }


          // then call the petsc function
          // that summarily allocates these
          // entries:
          ierr = MatMPIAIJSetPreallocationCSR (matrix,
                                               &rowstart_in_window[0],
                                               &colnums_in_window[0],
                                               0);
          AssertThrow (ierr == 0, ExcPETScError(ierr));

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
                ierr = MatSetValues (matrix, 1, &petsc_i,
                                     sparsity_pattern.row_length(i),
                                     &colnums_in_window[rowstart_in_window[i-local_row_start]],
                                     &values[0], INSERT_VALUES);
                AssertThrow (ierr == 0, ExcPETScError(ierr));
              }
          }

          compress (VectorOperation::insert);

          // set the dummy entries set above
          // back to zero
          *this = 0;
#endif // version <=2.3.3

          close_matrix (matrix);
          set_keep_zero_rows(matrix);
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
                                const DynamicSparsityPattern &,
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
                          const DynamicSparsityPattern &,
                          const std::vector<size_type> &,
                          const std::vector<size_type> &,
                          const unsigned int,
                          const bool);

    template void
    SparseMatrix::
    reinit (const IndexSet &,
            const IndexSet &,
            const SparsityPattern &,
            const MPI_Comm &);

    template void
    SparseMatrix::
    reinit (const IndexSet &,
            const IndexSet &,
            const DynamicSparsityPattern &,
            const MPI_Comm &);

    template void
    SparseMatrix::do_reinit (const SparsityPattern &,
                             const std::vector<size_type> &,
                             const std::vector<size_type> &,
                             const unsigned int ,
                             const bool);
    template void
    SparseMatrix::do_reinit (const DynamicSparsityPattern &,
                             const std::vector<size_type> &,
                             const std::vector<size_type> &,
                             const unsigned int ,
                             const bool);

    template void
    SparseMatrix::
    do_reinit (const IndexSet &,
               const IndexSet &,
               const SparsityPattern &);

    template void
    SparseMatrix::
    do_reinit (const IndexSet &,
               const IndexSet &,
               const DynamicSparsityPattern &);


    PetscScalar
    SparseMatrix::matrix_norm_square (const Vector &v) const
    {
      Vector tmp (v);
      vmult (tmp, v);
      // note, that v*tmp returns  sum_i conjugate(v)_i * tmp_i
      return v*tmp;
    }

    PetscScalar
    SparseMatrix::matrix_scalar_product (const Vector &u,
                                         const Vector &v) const
    {
      Vector tmp (v);
      vmult (tmp, v);
      // note, that v*tmp returns  sum_i conjugate(v)_i * tmp_i
      return u*tmp;
    }

    IndexSet
    SparseMatrix::locally_owned_domain_indices () const
    {
#if DEAL_II_PETSC_VERSION_LT(3,3,0)
      Assert(false,ExcNotImplemented());
      return IndexSet();
#else
      PetscInt n_rows, n_cols, n_loc_rows, n_loc_cols, min, max;
      PetscErrorCode ierr;

      ierr = MatGetSize (matrix, &n_rows, &n_cols);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      ierr = MatGetLocalSize(matrix, &n_loc_rows, &n_loc_cols);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      ierr = MatGetOwnershipRangeColumn(matrix, &min, &max);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      Assert(n_loc_cols==max-min, ExcMessage("PETSc is requiring non contiguous memory allocation."));

      IndexSet indices(n_cols);
      indices.add_range(min, max);
      indices.compress();

      return indices;
#endif
    }

    IndexSet
    SparseMatrix::locally_owned_range_indices () const
    {
#if DEAL_II_PETSC_VERSION_LT(3,3,0)
      Assert(false,ExcNotImplemented());
      return IndexSet();
#else
      PetscInt n_rows, n_cols, n_loc_rows, n_loc_cols, min, max;
      PetscErrorCode ierr;

      ierr = MatGetSize (matrix, &n_rows, &n_cols);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      ierr = MatGetLocalSize(matrix, &n_loc_rows, &n_loc_cols);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      ierr = MatGetOwnershipRange(matrix, &min, &max);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      Assert(n_loc_rows==max-min, ExcMessage("PETSc is requiring non contiguous memory allocation."));

      IndexSet indices(n_rows);
      indices.add_range(min, max);
      indices.compress();

      return indices;
#endif
    }

  }
}


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
