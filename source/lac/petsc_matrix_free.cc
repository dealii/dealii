// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/lac/petsc_matrix_free.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_compatibility.h>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  MatrixFree::MatrixFree()
  {
    const int m = 0;
    do_reinit(MPI_COMM_SELF, m, m, m, m);
  }



  MatrixFree::MatrixFree(const MPI_Comm     communicator,
                         const unsigned int m,
                         const unsigned int n,
                         const unsigned int local_rows,
                         const unsigned int local_columns)
  {
    do_reinit(communicator, m, n, local_rows, local_columns);
  }



  MatrixFree::MatrixFree(
    const MPI_Comm                   communicator,
    const unsigned int               m,
    const unsigned int               n,
    const std::vector<unsigned int> &local_rows_per_process,
    const std::vector<unsigned int> &local_columns_per_process,
    const unsigned int               this_process)
  {
    Assert(local_rows_per_process.size() == local_columns_per_process.size(),
           ExcDimensionMismatch(local_rows_per_process.size(),
                                local_columns_per_process.size()));
    Assert(this_process < local_rows_per_process.size(), ExcInternalError());

    do_reinit(communicator,
              m,
              n,
              local_rows_per_process[this_process],
              local_columns_per_process[this_process]);
  }



  MatrixFree::MatrixFree(const unsigned int m,
                         const unsigned int n,
                         const unsigned int local_rows,
                         const unsigned int local_columns)
  {
    do_reinit(MPI_COMM_WORLD, m, n, local_rows, local_columns);
  }



  MatrixFree::MatrixFree(
    const unsigned int               m,
    const unsigned int               n,
    const std::vector<unsigned int> &local_rows_per_process,
    const std::vector<unsigned int> &local_columns_per_process,
    const unsigned int               this_process)
  {
    Assert(local_rows_per_process.size() == local_columns_per_process.size(),
           ExcDimensionMismatch(local_rows_per_process.size(),
                                local_columns_per_process.size()));
    Assert(this_process < local_rows_per_process.size(), ExcInternalError());

    do_reinit(MPI_COMM_WORLD,
              m,
              n,
              local_rows_per_process[this_process],
              local_columns_per_process[this_process]);
  }



  void
  MatrixFree::reinit(const MPI_Comm     communicator,
                     const unsigned int m,
                     const unsigned int n,
                     const unsigned int local_rows,
                     const unsigned int local_columns)
  {
    // destroy the matrix and generate a new one
    const PetscErrorCode ierr = MatDestroy(&matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    do_reinit(communicator, m, n, local_rows, local_columns);
  }



  void
  MatrixFree::reinit(const MPI_Comm                   communicator,
                     const unsigned int               m,
                     const unsigned int               n,
                     const std::vector<unsigned int> &local_rows_per_process,
                     const std::vector<unsigned int> &local_columns_per_process,
                     const unsigned int               this_process)
  {
    Assert(local_rows_per_process.size() == local_columns_per_process.size(),
           ExcDimensionMismatch(local_rows_per_process.size(),
                                local_columns_per_process.size()));
    Assert(this_process < local_rows_per_process.size(), ExcInternalError());

    const PetscErrorCode ierr = MatDestroy(&matrix);
    AssertThrow(ierr != 0, ExcPETScError(ierr));

    do_reinit(communicator,
              m,
              n,
              local_rows_per_process[this_process],
              local_columns_per_process[this_process]);
  }



  void
  MatrixFree::reinit(const unsigned int m,
                     const unsigned int n,
                     const unsigned int local_rows,
                     const unsigned int local_columns)
  {
    reinit(this->get_mpi_communicator(), m, n, local_rows, local_columns);
  }



  void
  MatrixFree::reinit(const unsigned int               m,
                     const unsigned int               n,
                     const std::vector<unsigned int> &local_rows_per_process,
                     const std::vector<unsigned int> &local_columns_per_process,
                     const unsigned int               this_process)
  {
    reinit(this->get_mpi_communicator(),
           m,
           n,
           local_rows_per_process,
           local_columns_per_process,
           this_process);
  }



  void
  MatrixFree::clear()
  {
    const PetscErrorCode ierr = MatDestroy(&matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    const int m = 0;
    do_reinit(MPI_COMM_SELF, m, m, m, m);
  }



  void
  MatrixFree::vmult(Vec &dst, const Vec &src) const
  {
    // VectorBase permits us to manipulate, but not own, a Vec
    PETScWrappers::VectorBase x(src);
    PETScWrappers::VectorBase y(dst);

    // This is implemented by derived classes
    vmult(y, x);
  }



  int
  MatrixFree::matrix_free_mult(Mat A, Vec src, Vec dst)
  {
    // create a pointer to this MatrixFree
    // object and link the given matrix A
    // to the matrix-vector multiplication
    // of this MatrixFree object,
    void                *this_object;
    const PetscErrorCode ierr = MatShellGetContext(A, &this_object);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // call vmult of this object:
    reinterpret_cast<MatrixFree *>(this_object)->vmult(dst, src);

    return (0);
  }



  void
  MatrixFree::do_reinit(const MPI_Comm     communicator,
                        const unsigned int m,
                        const unsigned int n,
                        const unsigned int local_rows,
                        const unsigned int local_columns)
  {
    Assert(local_rows <= m, ExcDimensionMismatch(local_rows, m));
    Assert(local_columns <= n, ExcDimensionMismatch(local_columns, n));

    // create a PETSc MatShell matrix-type
    // object of dimension m x n and local size
    // local_rows x local_columns
    PetscErrorCode ierr = MatCreateShell(communicator,
                                         local_rows,
                                         local_columns,
                                         m,
                                         n,
                                         static_cast<void *>(this),
                                         &matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    // register the MatrixFree::matrix_free_mult function
    // as the matrix multiplication used by this matrix
    ierr = MatShellSetOperation(
      matrix,
      MATOP_MULT,
      reinterpret_cast<void (*)()>(
        &dealii::PETScWrappers::MatrixFree::matrix_free_mult));
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = MatSetFromOptions(matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }
} // namespace PETScWrappers


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
