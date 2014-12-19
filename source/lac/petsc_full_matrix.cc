// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2014 by the deal.II authors
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

#include <deal.II/lac/petsc_full_matrix.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/petsc_vector.h>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{

  FullMatrix::FullMatrix ()
  {
    // empty constructor generate an empty matrix
    do_reinit (0, 0);
  }

  FullMatrix::FullMatrix (const size_type m,
                          const size_type n)
  {
    do_reinit (m, n);
  }

  void
  FullMatrix::reinit (const size_type m,
                      const size_type n)
  {
    // get rid of old matrix and generate a
    // new one
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    const int ierr = MatDestroy (matrix);
#else
    const int ierr = MatDestroy (&matrix);
#endif
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    do_reinit (m, n);
  }

  void
  FullMatrix::do_reinit (const size_type m,
                         const size_type n)
  {
    // use the call sequence indicating only a maximal number of
    // elements per row for all rows globally
    const int ierr
      = MatCreateSeqDense (PETSC_COMM_SELF, m, n, PETSC_NULL,
                           &matrix);

    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }


  const MPI_Comm &
  FullMatrix::get_mpi_communicator () const
  {
    static const MPI_Comm communicator = MPI_COMM_SELF;
    return communicator;
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
