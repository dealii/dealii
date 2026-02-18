// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/lac/matrix_out.h>

#include <string>
#include <vector>


DEAL_II_NAMESPACE_OPEN


MatrixOut::Options::Options(const bool         show_absolute_values,
                            const unsigned int block_size,
                            const bool         discontinuous,
                            const bool         create_sparse_plot)
  : show_absolute_values(show_absolute_values)
  , block_size(block_size)
  , discontinuous(discontinuous)
  , create_sparse_plot(create_sparse_plot)
{}



const std::vector<MatrixOut::Patch> &
MatrixOut::get_patches() const
{
  return patches;
}



std::vector<std::string>
MatrixOut::get_dataset_names() const
{
  return std::vector<std::string>(1, name);
}

DEAL_II_NAMESPACE_CLOSE
