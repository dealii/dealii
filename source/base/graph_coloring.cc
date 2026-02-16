// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 1998 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/graph_coloring.h>

#include <deal.II/lac/sparsity_tools.h>

DEAL_II_NAMESPACE_OPEN

namespace GraphColoring
{
  unsigned int
  color_sparsity_pattern(const SparsityPattern     &sparsity_pattern,
                         std::vector<unsigned int> &color_indices)
  {
    return SparsityTools::color_sparsity_pattern(sparsity_pattern,
                                                 color_indices);
  }
} // namespace GraphColoring

DEAL_II_NAMESPACE_CLOSE
