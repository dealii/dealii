// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/distributed/cell_data_transfer.templates.h>

#ifdef DEAL_II_WITH_P4EST

DEAL_II_NAMESPACE_OPEN

// explicit instantiations
#  include "distributed/cell_data_transfer.inst"

DEAL_II_NAMESPACE_CLOSE

#endif /* DEAL_II_WITH_P4EST */
