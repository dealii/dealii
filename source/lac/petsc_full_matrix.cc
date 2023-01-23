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

#include <deal.II/lac/petsc_full_matrix.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_compatibility.h>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  FullMatrix::FullMatrix()
  {
    // empty constructor generate an empty matrix
    do_reinit(0, 0);
  }

  FullMatrix::FullMatrix(const size_type m, const size_type n)
  {
    do_reinit(m, n);
  }

  void
  FullMatrix::reinit(const size_type m, const size_type n)
  {
    // get rid of old matrix and generate a
    // new one
    const PetscErrorCode ierr = destroy_matrix(matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    do_reinit(m, n);
  }

  void
  FullMatrix::do_reinit(const size_type m, const size_type n)
  {
    // use the call sequence indicating only a maximal number of
    // elements per row for all rows globally
    const PetscErrorCode ierr =
      MatCreateSeqDense(PETSC_COMM_SELF, m, n, nullptr, &matrix);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }

} // namespace PETScWrappers


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
