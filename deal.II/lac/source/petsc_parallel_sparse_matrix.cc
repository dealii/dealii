//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/petsc_parallel_sparse_matrix.h>
#include <lac/petsc_vector.h>

#ifdef DEAL_II_USE_PETSC


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



    SparseMatrix::SparseMatrix (const MPI_Comm    &communicator,
                                const unsigned int m,
                                const unsigned int n,
                                const unsigned int local_rows,
                                const unsigned int local_columns,
                                const unsigned int n_nonzero_per_row,
                                const bool         is_symmetric)
                    :
                    communicator (communicator)
    {
      do_reinit (m, n, local_rows, local_columns,
                 n_nonzero_per_row, is_symmetric);
    }



    SparseMatrix::SparseMatrix (const MPI_Comm    &communicator,
                                const unsigned int m,
                                const unsigned int n,
                                const unsigned int local_rows,
                                const unsigned int local_columns,
                                const std::vector<unsigned int> &row_lengths,
                                const bool         is_symmetric)
                    :
                    communicator (communicator)
    {
      do_reinit (m, n, local_rows, local_columns,
                 row_lengths, is_symmetric);
    }



    SparseMatrix::
    SparseMatrix (const MPI_Comm                  &communicator,
                  const CompressedSparsityPattern &sparsity_pattern,
                  const std::vector<unsigned int> &local_rows_per_process,
                  const std::vector<unsigned int> &local_columns_per_process,
                  const unsigned int               this_process,
                  const bool                       preset_nonzero_locations)
                    :
                    communicator (communicator)
    {
      do_reinit (sparsity_pattern, local_rows_per_process,
                 local_columns_per_process, this_process,
                 preset_nonzero_locations);
    }
    


    SparseMatrix &
    SparseMatrix::operator = (const double d)
    {
      MatrixBase::operator = (d);
      return *this;
    }

    

    void
    SparseMatrix::reinit (const MPI_Comm    &communicator,
                          const unsigned int m,
                          const unsigned int n,
                          const unsigned int local_rows,
                          const unsigned int local_columns,
                          const unsigned int n_nonzero_per_row,
                          const bool         is_symmetric)
    {
      this->communicator = communicator;

                                       // get rid of old matrix and generate a
                                       // new one
      const int ierr = MatDestroy (matrix);
      AssertThrow (ierr == 0, ExcPETScError(ierr));    
    
      do_reinit (m, n, local_rows, local_columns,
                 n_nonzero_per_row, is_symmetric);
    }



    void
    SparseMatrix::reinit (const MPI_Comm    &communicator,
                          const unsigned int m,
                          const unsigned int n,
                          const unsigned int local_rows,
                          const unsigned int local_columns,
                          const std::vector<unsigned int> &row_lengths,
                          const bool         is_symmetric)
    {
      this->communicator = communicator;
      
                                       // get rid of old matrix and generate a
                                       // new one
      const int ierr = MatDestroy (matrix);
      AssertThrow (ierr == 0, ExcPETScError(ierr));    

      do_reinit (m, n, local_rows, local_columns, row_lengths, is_symmetric);
    }  



    void
    SparseMatrix::
    reinit (const MPI_Comm                  &communicator,
            const CompressedSparsityPattern &sparsity_pattern,
            const std::vector<unsigned int> &local_rows_per_process,
            const std::vector<unsigned int> &local_columns_per_process,
            const unsigned int               this_process,
            const bool                       preset_nonzero_locations)
    {
      this->communicator = communicator;
      
                                       // get rid of old matrix and generate a
                                       // new one
      const int ierr = MatDestroy (matrix);
      AssertThrow (ierr == 0, ExcPETScError(ierr));    

      do_reinit (sparsity_pattern, local_rows_per_process,
                 local_columns_per_process, this_process,
                 preset_nonzero_locations);
    }

    

    void
    SparseMatrix::do_reinit (const unsigned int m,
                             const unsigned int n,
                             const unsigned int local_rows,
                             const unsigned int local_columns,
                             const unsigned int n_nonzero_per_row,
                             const bool         is_symmetric)
    {
      Assert (local_rows <= m, ExcLocalRowsTooLarge (local_rows, m));
      
                                       // use the call sequence indicating only
                                       // a maximal number of elements per row
                                       // for all rows globally
      const int ierr
        = MatCreateMPIAIJ(communicator,
			  local_rows, local_columns,
			  m, n,
                          n_nonzero_per_row, 0, 0, 0,
                          &matrix);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

                                       // set symmetric flag, if so requested
      if (is_symmetric == true)
        {
          const int ierr
            = MatSetOption (matrix, MAT_SYMMETRIC);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
        }
    }



    void
    SparseMatrix::do_reinit (const unsigned int m,
                             const unsigned int n,
                             const unsigned int local_rows,
                             const unsigned int local_columns,
                             const std::vector<unsigned int> &row_lengths,
                             const bool         is_symmetric)
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
      for (unsigned int i=0; i<row_lengths.size(); ++i)
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
      const std::vector<signed int> int_row_lengths (row_lengths.begin(),
                                                     row_lengths.end());

//TODO: There must be a significantly better way to provide information about the off-diagonal blocks of the matrix. this way, petsc keeps allocating tiny chunks of memory, and gets completely hung up over this
      const int ierr
        = MatCreateMPIAIJ(communicator,
			  local_rows, local_columns,
			  m, n,
                          0, &int_row_lengths[0], 0, 0,
                          &matrix);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

                                       // set symmetric flag, if so requested
      if (is_symmetric == true)
        {
          const int ierr
            = MatSetOption (matrix, MAT_SYMMETRIC);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
        }    
    }



    void
    SparseMatrix::
    do_reinit (const CompressedSparsityPattern &sparsity_pattern,
               const std::vector<unsigned int> &local_rows_per_process,
               const std::vector<unsigned int> &local_columns_per_process,
               const unsigned int               this_process,
               const bool                       preset_nonzero_locations)
    {
      Assert (local_rows_per_process.size() == local_columns_per_process.size(),
              ExcInvalidSizes (local_rows_per_process.size(),
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
      unsigned int local_row_start = 0;
      unsigned int local_col_start = 0;
      for (unsigned int p=0; p<this_process; ++p)
        {
          local_row_start += local_rows_per_process[p];
          local_col_start += local_columns_per_process[p];
        }
      const unsigned int
        local_row_end = local_row_start + local_rows_per_process[this_process];
      const unsigned int
        local_col_end = local_col_start + local_columns_per_process[this_process];

                                       // then count the elements in- and
                                       // out-of-window for the rows we own
      std::vector<int> row_lengths_in_window (local_row_end - local_row_start);
      std::vector<int> row_lengths_out_of_window (local_row_end - local_row_start);
      for (unsigned int row = local_row_start; row<local_row_end; ++row)
        for (unsigned int c=0; c<sparsity_pattern.row_length(row); ++c)
          {
            const unsigned int column = sparsity_pattern.column_number(row,c);
            
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
          std::vector<int> row_entries;
          std::vector<PetscScalar> row_values;
          for (unsigned int i=0; i<sparsity_pattern.n_rows(); ++i)
            {
              row_entries.resize (sparsity_pattern.row_length(i));
              row_values.resize (sparsity_pattern.row_length(i), 0.0);
              for (unsigned int j=0; j<sparsity_pattern.row_length(i); ++j)
                row_entries[j] = sparsity_pattern.column_number (i,j);
              
              const int int_row = i;
              MatSetValues (matrix, 1, &int_row,
                            sparsity_pattern.row_length(i), &row_entries[0],
                            &row_values[0], INSERT_VALUES);
            }
          compress ();

//TODO: We should use the appropriate option with MatSetOption here indicating that we do not intend to add any further matrix entries. Maybe this would accelerate some things as well
        }
    }
    
  }
}


#else
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // DEAL_II_USE_PETSC
