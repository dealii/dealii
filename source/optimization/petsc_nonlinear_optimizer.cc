// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef DOXYGEN

#  include <deal.II/base/config.h>

#  if defined(DEAL_II_WITH_PETSC) && DEAL_II_PETSC_VERSION_GTE(3, 13, 0)

#    include <deal.II/lac/petsc_block_sparse_matrix.h>

#    include <deal.II/optimization/petsc_nonlinear_optimizer.templates.h>


#  endif
#endif

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
#  if defined(DEAL_II_WITH_PETSC) && DEAL_II_PETSC_VERSION_GTE(3, 13, 0)

template class PETScWrappers::NonlinearOptimizer<>;
template class PETScWrappers::NonlinearOptimizer<
  PETScWrappers::MPI::Vector,
  PETScWrappers::MPI::SparseMatrix>;
template class PETScWrappers::NonlinearOptimizer<
  PETScWrappers::MPI::BlockVector,
  PETScWrappers::MPI::BlockSparseMatrix>;



#  endif
#endif

DEAL_II_NAMESPACE_CLOSE
