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


#include <deal.II/lac/petsc_full_matrix.h>

#ifdef DEAL_II_USE_PETSC

#  include <deal.II/lac/petsc_vector.h>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  FullMatrix::FullMatrix (const unsigned int m,
                          const unsigned int n)
  {
    const int ierr
      = MatCreateSeqDense(PETSC_COMM_SELF, m, n, PETSC_NULL,
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

#endif // DEAL_II_USE_PETSC
