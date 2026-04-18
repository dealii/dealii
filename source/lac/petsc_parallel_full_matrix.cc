// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2004 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/lac/petsc_full_matrix.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/base/mpi.h>

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_compatibility.h>
#  include <deal.II/lac/petsc_vector.h>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  namespace MPI
  {
    FullMatrix::FullMatrix()
    {
      // just like for vectors and sparse matrices: since we
      // create an empty matrix, we can as
      // well make it sequential
      const int            m = 0, n = 0;
      const PetscErrorCode ierr =
        MatCreateSeqDense(PETSC_COMM_SELF, m, n, nullptr, &matrix);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
    }



    FullMatrix::FullMatrix(const Mat &A)
      : MatrixBase(A)
    {}



    FullMatrix::~FullMatrix()
    {
      PetscErrorCode ierr = MatDestroy(&matrix);
      AssertNothrow(ierr == 0, ExcPETScError(ierr));
    }



    FullMatrix::FullMatrix(
      const MPI_Comm                communicator,
      const size_type               m,
      const size_type               n,
      const std::vector<size_type> &local_rows_per_process,
      const std::vector<size_type> &local_columns_per_process,
      const unsigned int            this_process)
    {
      do_reinit(communicator,
                m,
                n,
                local_rows_per_process,
                local_columns_per_process,
                this_process);
    }



    void
    FullMatrix::reinit(const MPI_Comm                communicator,
                       const size_type               m,
                       const size_type               n,
                       const std::vector<size_type> &local_rows_per_process,
                       const std::vector<size_type> &local_columns_per_process,
                       const unsigned int            this_process)
    {
      // get rid of old matrix and generate a new one
      const PetscErrorCode ierr = MatDestroy(&matrix);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      do_reinit(communicator,
                m,
                n,
                local_rows_per_process,
                local_columns_per_process,
                this_process);
    }



    void
    FullMatrix::do_reinit(
      const MPI_Comm                communicator,
      const size_type               m,
      const size_type               n,
      const std::vector<size_type> &local_rows_per_process,
      const std::vector<size_type> &local_columns_per_process,
      const unsigned int            this_process)
    {
      Assert(this_process < local_rows_per_process.size(), ExcInternalError());

      // convert dimensions into PETScInt
      AssertThrowIntegerConversion(static_cast<PetscInt>(m), m);
      AssertThrowIntegerConversion(static_cast<PetscInt>(n), n);

      // create a matrix
      // this will create by default a sequential dense matrix(MATSEQDENSE) if
      // constructed with a single process communicator and a distributed dense
      // matrix (MATMPIDENSE) if constructed with more than one.
      PetscErrorCode ierr =
        MatCreateDense(communicator,
                       local_rows_per_process[this_process],
                       local_columns_per_process[this_process],
                       m,
                       n,
                       nullptr,
                       &matrix);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
    }



  } // namespace MPI

} // namespace PETScWrappers
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
