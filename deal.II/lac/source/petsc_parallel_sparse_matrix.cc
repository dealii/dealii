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
                  const unsigned int               local_rows,
                  const unsigned int               local_columns,
                  const bool                       preset_nonzero_locations)
                    :
                    communicator (communicator)
    {
      do_reinit (sparsity_pattern, local_rows, local_columns,
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
            const unsigned int               local_rows,
            const unsigned int               local_columns,
            const bool                       preset_nonzero_locations)
    {
      this->communicator = communicator;
      
                                       // get rid of old matrix and generate a
                                       // new one
      const int ierr = MatDestroy (matrix);
      AssertThrow (ierr == 0, ExcPETScError(ierr));    

      do_reinit (sparsity_pattern, local_rows, local_columns,
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
//TODO[WB]: We should do better by providing ways to tell PETSc how to partition the columns
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
                                       
//TODO[WB]: We should do better by providing ways to tell PETSc how to partition the columns
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
    SparseMatrix::do_reinit (const CompressedSparsityPattern &sparsity_pattern,
                             const unsigned int               local_rows,
                             const unsigned int               local_columns,
                             const bool preset_nonzero_locations)
    {
      std::vector<unsigned int> row_lengths (sparsity_pattern.n_rows());
      for (unsigned int i=0; i<sparsity_pattern.n_rows(); ++i)
        row_lengths[i] = sparsity_pattern.row_length (i);

      do_reinit (sparsity_pattern.n_rows(),
                 sparsity_pattern.n_cols(),
                 local_rows,
                 local_columns,
                 row_lengths, false);

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
              row_entries.resize (row_lengths[i]);
              row_values.resize (row_lengths[i], 0.0);
              for (unsigned int j=0; j<row_lengths[i]; ++j)
                row_entries[j] = sparsity_pattern.column_number (i,j);
              
              const int int_row = i;
              MatSetValues (matrix, 1, &int_row,
                            row_lengths[i], &row_entries[0],
                            &row_values[0], INSERT_VALUES);
            }
          compress ();
        }
    }
    
  }
}


#else
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // DEAL_II_USE_PETSC
