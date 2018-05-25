//-----------------------------------------------------------
//
//    Copyright (C) 2017 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE at
//    the top level of the deal.II distribution.
//
//-----------------------------------------------------------

#ifndef dealii_sundials_copy_h
#define dealii_sundials_copy_h

#include <deal.II/base/config.h>
#ifdef DEAL_II_WITH_SUNDIALS

#  include <sundials/sundials_nvector.h>
#  ifdef DEAL_II_WITH_MPI
#    include <nvector/nvector_parallel.h>
#  endif
#  include <deal.II/lac/block_vector.h>
#  include <deal.II/lac/vector.h>

#  include <nvector/nvector_serial.h>

#  ifdef DEAL_II_WITH_TRILINOS
#    include <deal.II/lac/trilinos_parallel_block_vector.h>
#    include <deal.II/lac/trilinos_vector.h>
#  endif

#  ifdef DEAL_II_WITH_PETSC
#    include <deal.II/lac/petsc_parallel_block_vector.h>
#    include <deal.II/lac/petsc_parallel_vector.h>
#  endif

DEAL_II_NAMESPACE_OPEN
namespace SUNDIALS
{
  namespace internal
  {
    // The following internal functions are used by SUNDIALS wrappers to copy
    // to and from deal.II vector types.
#  ifdef DEAL_II_WITH_MPI

#    ifdef DEAL_II_WITH_TRILINOS
    void
    copy(TrilinosWrappers::MPI::Vector &dst, const N_Vector &src);
    void
    copy(N_Vector &dst, const TrilinosWrappers::MPI::Vector &src);
    void
    copy(TrilinosWrappers::MPI::BlockVector &dst, const N_Vector &src);
    void
    copy(N_Vector &dst, const TrilinosWrappers::MPI::BlockVector &src);
#    endif // DEAL_II_WITH_TRILINOS

#    ifdef DEAL_II_WITH_PETSC
#      ifndef PETSC_USE_COMPLEX
    void
    copy(PETScWrappers::MPI::Vector &dst, const N_Vector &src);
    void
    copy(N_Vector &dst, const PETScWrappers::MPI::Vector &src);
    void
    copy(PETScWrappers::MPI::BlockVector &dst, const N_Vector &src);
    void
    copy(N_Vector &dst, const PETScWrappers::MPI::BlockVector &src);
#      endif // PETSC_USE_COMPLEX
#    endif   // DEAL_II_WITH_PETSC

#  endif

    void
    copy(BlockVector<double> &dst, const N_Vector &src);
    void
    copy(N_Vector &dst, const BlockVector<double> &src);

    void
    copy(Vector<double> &dst, const N_Vector &src);
    void
    copy(N_Vector &dst, const Vector<double> &src);
  } // namespace internal
} // namespace SUNDIALS
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SUNDIALS
#endif // dealii_sundials_copy_h
