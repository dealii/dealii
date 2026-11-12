// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2024 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/lac/trilinos_tpetra_block_sparse_matrix.templates.h>


#endif // DEAL_II_TRILINOS_WITH_TPETRA

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  ifndef DOXYGEN

#    include "lac/trilinos_tpetra_block_sparse_matrix.inst"

#  endif // DOXYGEN


#endif // DEAL_II_TRILINOS_WITH_TPETRA
DEAL_II_NAMESPACE_CLOSE
