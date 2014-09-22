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

#include <deal.II/lac/petsc_sparse_matrix.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/sparsity_pattern.h>
#  include <deal.II/lac/compressed_sparsity_pattern.h>
#  include <deal.II/lac/compressed_simple_sparsity_pattern.h>

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



  SparseMatrix::SparseMatrix (const size_type m,
                              const size_type n,
                              const size_type n_nonzero_per_row,
                              const bool      is_symmetric)
  {
    do_reinit (m, n, n_nonzero_per_row, is_symmetric);
  }



  SparseMatrix::SparseMatrix (const size_type               m,
                              const size_type               n,
                              const std::vector<size_type> &row_lengths,
                              const bool                    is_symmetric)
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
  SparseMatrix::reinit (const size_type m,
                        const size_type n,
                        const size_type n_nonzero_per_row,
                        const bool      is_symmetric)
  {
    // get rid of old matrix and generate a
    // new one
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    const int ierr = MatDestroy (matrix);
#else
    const int ierr = MatDestroy (&matrix);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    do_reinit (m, n, n_nonzero_per_row, is_symmetric);
  }



  void
  SparseMatrix::reinit (const size_type               m,
                        const size_type               n,
                        const std::vector<size_type> &row_lengths,
                        const bool                    is_symmetric)
  {
    // get rid of old matrix and generate a
    // new one
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    const int ierr = MatDestroy (matrix);
#else
    const int ierr = MatDestroy (&matrix);
#endif
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
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    const int ierr = MatDestroy (matrix);
#else
    const int ierr = MatDestroy (&matrix);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    do_reinit (sparsity_pattern, preset_nonzero_locations);
  }



  const MPI_Comm &
  SparseMatrix::get_mpi_communicator () const
  {
    static MPI_Comm comm;
    PetscObjectGetComm((PetscObject)matrix, &comm);
    return comm;
  }



  void
  SparseMatrix::do_reinit (const size_type m,
                           const size_type n,
                           const size_type n_nonzero_per_row,
                           const bool      is_symmetric)
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
  SparseMatrix::do_reinit (const size_type               m,
                           const size_type               n,
                           const std::vector<size_type> &row_lengths,
                           const bool                    is_symmetric)
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
    const std::vector<PetscInt>
    int_row_lengths (row_lengths.begin(), row_lengths.end());

    const int ierr
      = MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, 0,
                        &int_row_lengths[0], &matrix);
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
  SparseMatrix::do_reinit (const SparsityType &sparsity_pattern,
                           const bool          preset_nonzero_locations)
  {
    std::vector<size_type> row_lengths (sparsity_pattern.n_rows());
    for (size_type i=0; i<sparsity_pattern.n_rows(); ++i)
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
        std::vector<PetscInt>    row_entries;
        std::vector<PetscScalar> row_values;
        for (size_type i=0; i<sparsity_pattern.n_rows(); ++i)
          {
            row_entries.resize (row_lengths[i]);
            row_values.resize (row_lengths[i], 0.0);
            for (size_type j=0; j<row_lengths[i]; ++j)
              row_entries[j] = sparsity_pattern.column_number (i,j);

            const PetscInt int_row = i;
            MatSetValues (matrix, 1, &int_row,
                          row_lengths[i], &row_entries[0],
                          &row_values[0], INSERT_VALUES);
          }
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

  PetscScalar
  SparseMatrix::matrix_norm_square (const VectorBase &v) const
  {
    Vector tmp (v.size());
    vmult (tmp, v);
    return tmp*v;
  }

  PetscScalar
  SparseMatrix::matrix_scalar_product (const VectorBase &u,
                                       const VectorBase &v) const
  {
    Vector tmp (v.size());
    vmult (tmp, v);
    return u*tmp;
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
