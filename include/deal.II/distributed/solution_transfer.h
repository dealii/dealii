// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_distributed_solution_transfer_h
#define dealii_distributed_solution_transfer_h

#include <deal.II/base/config.h>

#include <deal.II/numerics/solution_transfer.h>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    /**
     * @deprecated Use dealii::SolutionTransfer instead.
     */
    template <int dim, typename VectorType, int spacedim = dim>
    using SolutionTransfer DEAL_II_DEPRECATED =
      dealii::SolutionTransfer<dim, VectorType, spacedim>;
  } // namespace distributed
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif
