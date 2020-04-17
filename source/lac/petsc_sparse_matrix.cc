// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
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

#include <deal.II/lac/petsc_sparse_matrix.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/dynamic_sparsity_pattern.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_compatibility.h>
#  include <deal.II/lac/petsc_vector_base.h>
#  include <deal.II/lac/sparsity_pattern.h>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  SparseMatrix::SparseMatrix()
  {
    const int            m = 0, n = 0, n_nonzero_per_row = 0;
    const PetscErrorCode ierr = MatCreateSeqAIJ(
      PETSC_COMM_SELF, m, n, n_nonzero_per_row, nullptr, &matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  SparseMatrix::SparseMatrix(const size_type m,
                             const size_type n,
                             const size_type n_nonzero_per_row,
                             const bool      is_symmetric)
  {
    do_reinit(m, n, n_nonzero_per_row, is_symmetric);
  }



  SparseMatrix::SparseMatrix(const size_type               m,
                             const size_type               n,
                             const std::vector<size_type> &row_lengths,
                             const bool                    is_symmetric)
  {
    do_reinit(m, n, row_lengths, is_symmetric);
  }



  template <typename SparsityPatternType>
  SparseMatrix::SparseMatrix(const SparsityPatternType &sparsity_pattern,
                             const bool preset_nonzero_locations)
  {
    do_reinit(sparsity_pattern, preset_nonzero_locations);
  }



  SparseMatrix &
  SparseMatrix::operator=(const double d)
  {
    MatrixBase::operator=(d);
    return *this;
  }



  void
  SparseMatrix::reinit(const size_type m,
                       const size_type n,
                       const size_type n_nonzero_per_row,
                       const bool      is_symmetric)
  {
    // get rid of old matrix and generate a
    // new one
    const PetscErrorCode ierr = destroy_matrix(matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    do_reinit(m, n, n_nonzero_per_row, is_symmetric);
  }



  void
  SparseMatrix::reinit(const size_type               m,
                       const size_type               n,
                       const std::vector<size_type> &row_lengths,
                       const bool                    is_symmetric)
  {
    // get rid of old matrix and generate a
    // new one
    const PetscErrorCode ierr = destroy_matrix(matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    do_reinit(m, n, row_lengths, is_symmetric);
  }



  template <typename SparsityPatternType>
  void
  SparseMatrix::reinit(const SparsityPatternType &sparsity_pattern,
                       const bool                 preset_nonzero_locations)
  {
    // get rid of old matrix and generate a
    // new one
    const PetscErrorCode ierr = destroy_matrix(matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    do_reinit(sparsity_pattern, preset_nonzero_locations);
  }



  const MPI_Comm &
  SparseMatrix::get_mpi_communicator() const
  {
    static MPI_Comm      comm;
    const PetscErrorCode ierr =
      PetscObjectGetComm(reinterpret_cast<PetscObject>(matrix), &comm);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    return comm;
  }



  void
  SparseMatrix::do_reinit(const size_type m,
                          const size_type n,
                          const size_type n_nonzero_per_row,
                          const bool      is_symmetric)
  {
    // use the call sequence indicating only
    // a maximal number of elements per row
    // for all rows globally
    const PetscErrorCode ierr = MatCreateSeqAIJ(
      PETSC_COMM_SELF, m, n, n_nonzero_per_row, nullptr, &matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // set symmetric flag, if so requested
    if (is_symmetric == true)
      {
        set_matrix_option(matrix, MAT_SYMMETRIC, PETSC_TRUE);
      }
  }



  void
  SparseMatrix::do_reinit(const size_type               m,
                          const size_type               n,
                          const std::vector<size_type> &row_lengths,
                          const bool                    is_symmetric)
  {
    Assert(row_lengths.size() == m,
           ExcDimensionMismatch(row_lengths.size(), m));

    // use the call sequence indicating a
    // maximal number of elements for each
    // row individually. annoyingly, we
    // always use unsigned ints for cases
    // like this, while PETSc wants to see
    // signed integers. so we have to
    // convert, unless we want to play dirty
    // tricks with conversions of pointers
    const std::vector<PetscInt> int_row_lengths(row_lengths.begin(),
                                                row_lengths.end());

    const PetscErrorCode ierr = MatCreateSeqAIJ(
      PETSC_COMM_SELF, m, n, 0, int_row_lengths.data(), &matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // set symmetric flag, if so requested
    if (is_symmetric == true)
      {
        set_matrix_option(matrix, MAT_SYMMETRIC, PETSC_TRUE);
      }
  }



  template <typename SparsityPatternType>
  void
  SparseMatrix::do_reinit(const SparsityPatternType &sparsity_pattern,
                          const bool                 preset_nonzero_locations)
  {
    std::vector<size_type> row_lengths(sparsity_pattern.n_rows());
    for (size_type i = 0; i < sparsity_pattern.n_rows(); ++i)
      row_lengths[i] = sparsity_pattern.row_length(i);

    do_reinit(sparsity_pattern.n_rows(),
              sparsity_pattern.n_cols(),
              row_lengths,
              false);

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
        for (size_type i = 0; i < sparsity_pattern.n_rows(); ++i)
          {
            row_entries.resize(row_lengths[i]);
            row_values.resize(row_lengths[i], 0.0);
            for (size_type j = 0; j < row_lengths[i]; ++j)
              row_entries[j] = sparsity_pattern.column_number(i, j);

            const PetscInt       int_row = i;
            const PetscErrorCode ierr    = MatSetValues(matrix,
                                                     1,
                                                     &int_row,
                                                     row_lengths[i],
                                                     row_entries.data(),
                                                     row_values.data(),
                                                     INSERT_VALUES);
            AssertThrow(ierr == 0, ExcPETScError(ierr));
          }
        compress(VectorOperation::insert);

        close_matrix(matrix);
        set_keep_zero_rows(matrix);
      }
  }

  size_t
  SparseMatrix::m() const
  {
    PetscInt             m, n;
    const PetscErrorCode ierr = MatGetSize(matrix, &m, &n);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return m;
  }

  size_t
  SparseMatrix::n() const
  {
    PetscInt             m, n;
    const PetscErrorCode ierr = MatGetSize(matrix, &m, &n);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return n;
  }

  void
  SparseMatrix::mmult(SparseMatrix &      C,
                      const SparseMatrix &B,
                      const MPI::Vector & V) const
  {
    // Simply forward to the protected member function of the base class
    // that takes abstract matrix and vector arguments (to which the compiler
    // automatically casts the arguments).
    MatrixBase::mmult(C, B, V);
  }

  void
  SparseMatrix::Tmmult(SparseMatrix &      C,
                       const SparseMatrix &B,
                       const MPI::Vector & V) const
  {
    // Simply forward to the protected member function of the base class
    // that takes abstract matrix and vector arguments (to which the compiler
    // automatically casts the arguments).
    MatrixBase::Tmmult(C, B, V);
  }

#  ifndef DOXYGEN
  // Explicit instantiations
  //
  template SparseMatrix::SparseMatrix(const SparsityPattern &, const bool);
  template SparseMatrix::SparseMatrix(const DynamicSparsityPattern &,
                                      const bool);

  template void
  SparseMatrix::reinit(const SparsityPattern &, const bool);
  template void
  SparseMatrix::reinit(const DynamicSparsityPattern &, const bool);

  template void
  SparseMatrix::do_reinit(const SparsityPattern &, const bool);
  template void
  SparseMatrix::do_reinit(const DynamicSparsityPattern &, const bool);
#  endif
} // namespace PETScWrappers


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
