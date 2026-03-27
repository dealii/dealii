// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/lac/trilinos_tpetra_sparse_matrix.templates.h>

DEAL_II_NAMESPACE_OPEN

#  ifndef DOXYGEN

#    include "lac/trilinos_tpetra_sparse_matrix.inst"

namespace LinearAlgebra::TpetraWrappers
{
  template void
  SparseMatrix<double, MemorySpace::Host>::set(const size_type,
                                               const size_type,
                                               const size_type *,
                                               const float *,
                                               bool);

  template void
  SparseMatrix<double, MemorySpace::Default>::set(const size_type,
                                                  const size_type,
                                                  const size_type *,
                                                  const float *,
                                                  bool);
} // namespace LinearAlgebra::TpetraWrappers

#  endif

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA
