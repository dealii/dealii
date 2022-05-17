// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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

#ifndef DOXYGEN

#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_PETSC

#    include <deal.II/lac/petsc_block_sparse_matrix.h>
#    include <deal.II/lac/petsc_ts.templates.h>

DEAL_II_NAMESPACE_OPEN


template class PETScWrappers::TimeStepper<>;
template class PETScWrappers::TimeStepper<PETScWrappers::MPI::Vector>;
template class PETScWrappers::TimeStepper<PETScWrappers::MPI::BlockVector>;
template class PETScWrappers::TimeStepper<PETScWrappers::MPI::Vector,
                                          PETScWrappers::MPI::SparseMatrix>;
template class PETScWrappers::TimeStepper<
  PETScWrappers::MPI::BlockVector,
  PETScWrappers::MPI::BlockSparseMatrix>;


DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_PETSC
#endif
