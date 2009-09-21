//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2006, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/petsc_sparse_matrix.h>
#include <lac/petsc_vector.h>

#include <lac/sparsity_pattern.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/compressed_simple_sparsity_pattern.h>

#ifdef DEAL_II_USE_PETSC

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  SparseMatrix::SparseMatrix ()
  {
    const int m=0, n=0, n_nonzero_per_row=0;
    const int ierr
      = MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, n_nonzero_per_row,
                        0, &matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  SparseMatrix::SparseMatrix (const unsigned int m,
                              const unsigned int n,
                              const unsigned int n_nonzero_per_row,
                              const bool         is_symmetric)
  {
    do_reinit (m, n, n_nonzero_per_row, is_symmetric);
  }



  SparseMatrix::SparseMatrix (const unsigned int m,
                              const unsigned int n,
                              const std::vector<unsigned int> &row_lengths,
                              const bool         is_symmetric)
  {
    do_reinit (m, n, row_lengths, is_symmetric);
  }



  template <typename SparsityType>
  SparseMatrix::
  SparseMatrix (const SparsityType &sparsity_pattern,
                const bool          preset_nonzero_locations)
  {
    do_reinit (sparsity_pattern, preset_nonzero_locations);
  }



  SparseMatrix &
  SparseMatrix::operator = (const double d)
  {
    MatrixBase::operator = (d);
    return *this;
  }



  void
  SparseMatrix::reinit (const unsigned int m,
                        const unsigned int n,
                        const unsigned int n_nonzero_per_row,
                        const bool         is_symmetric)
  {
                                     // get rid of old matrix and generate a
                                     // new one
    const int ierr = MatDestroy (matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    do_reinit (m, n, n_nonzero_per_row, is_symmetric);
  }



  void
  SparseMatrix::reinit (const unsigned int m,
                        const unsigned int n,
                        const std::vector<unsigned int> &row_lengths,
                        const bool         is_symmetric)
  {
                                     // get rid of old matrix and generate a
                                     // new one
    const int ierr = MatDestroy (matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    do_reinit (m, n, row_lengths, is_symmetric);
  }



  template <typename SparsityType>
  void
  SparseMatrix::
  reinit (const SparsityType &sparsity_pattern,
          const bool          preset_nonzero_locations)
  {
                                     // get rid of old matrix and generate a
                                     // new one
    const int ierr = MatDestroy (matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    do_reinit (sparsity_pattern, preset_nonzero_locations);
  }



  const MPI_Comm &
  SparseMatrix::get_mpi_communicator () const
  {
    static const MPI_Comm communicator = MPI_COMM_SELF;
    return communicator;
  }



  void
  SparseMatrix::do_reinit (const unsigned int m,
                           const unsigned int n,
                           const unsigned int n_nonzero_per_row,
                           const bool         is_symmetric)
  {
                                     // use the call sequence indicating only
                                     // a maximal number of elements per row
                                     // for all rows globally
    const int ierr
      = MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, n_nonzero_per_row,
                        0, &matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

                                     // set symmetric flag, if so requested
    if (is_symmetric == true)
      {
#if (PETSC_VERSION_MAJOR <= 2)
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
  SparseMatrix::do_reinit (const unsigned int m,
                           const unsigned int n,
                           const std::vector<unsigned int> &row_lengths,
                           const bool         is_symmetric)
  {
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
#ifdef PETSC_USE_64BIT_INDICES
    const std::vector<PetscInt>
#else
    const std::vector<int>
#endif
      int_row_lengths (row_lengths.begin(), row_lengths.end());
    const int ierr
      = MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, 0,
                        &int_row_lengths[0], &matrix);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

                                     // set symmetric flag, if so requested
    if (is_symmetric == true)
      {
#if (PETSC_VERSION_MAJOR <= 2)
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
  SparseMatrix::do_reinit (const SparsityType &sparsity_pattern,
                           const bool          preset_nonzero_locations)
  {
    std::vector<unsigned int> row_lengths (sparsity_pattern.n_rows());
    for (unsigned int i=0; i<sparsity_pattern.n_rows(); ++i)
      row_lengths[i] = sparsity_pattern.row_length (i);

    do_reinit (sparsity_pattern.n_rows(),
               sparsity_pattern.n_cols(),
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
#ifdef PETSC_USE_64BIT_INDICES
	std::vector<PetscInt>
#else
	std::vector<int>
#endif
	  row_entries;
        std::vector<PetscScalar> row_values;
        for (unsigned int i=0; i<sparsity_pattern.n_rows(); ++i)
          {
            row_entries.resize (row_lengths[i]);
            row_values.resize (row_lengths[i], 0.0);
            for (unsigned int j=0; j<row_lengths[i]; ++j)
              row_entries[j] = sparsity_pattern.column_number (i,j);

#ifdef PETSC_USE_64BIT_INDICES
	    const PetscInt
#else
	    const int
#endif
	      int_row = i;
            MatSetValues (matrix, 1, &int_row,
                          row_lengths[i], &row_entries[0],
                          &row_values[0], INSERT_VALUES);
          }
        compress ();
      }

				        // In the end, tell the matrix that
				        // it should not expect any new
				        // entries.
#if (PETSC_VERSION_MAJOR <= 2)
    const int ierr =
      MatSetOption (matrix, MAT_NO_NEW_NONZERO_LOCATIONS);
#else
    const int ierr =
      MatSetOption (matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }


  // Explicit instantiations
  //
  template
  SparseMatrix::SparseMatrix (const SparsityPattern &,
			      const bool);
  template
  SparseMatrix::SparseMatrix (const CompressedSparsityPattern &,
			      const bool);
  template
  SparseMatrix::SparseMatrix (const CompressedSimpleSparsityPattern &,
			      const bool);

  template void
  SparseMatrix::reinit (const SparsityPattern &,
			const bool);
  template void
  SparseMatrix::reinit (const CompressedSparsityPattern &,
			const bool);
  template void
  SparseMatrix::reinit (const CompressedSimpleSparsityPattern &,
			const bool);

  template void
  SparseMatrix::do_reinit (const SparsityPattern &,
                           const bool);
  template void
  SparseMatrix::do_reinit (const CompressedSparsityPattern &,
                           const bool);
  template void
  SparseMatrix::do_reinit (const CompressedSimpleSparsityPattern &,
                           const bool);
}


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_PETSC
