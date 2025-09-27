// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
    const PetscErrorCode ierr = MatDestroy(&matrix);
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
