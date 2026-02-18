// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2003 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/lac/block_sparse_matrix_ez.h>
#include <deal.II/lac/block_sparse_matrix_ez.templates.h>

DEAL_II_NAMESPACE_OPEN

// explicit instantiations
template class BlockSparseMatrixEZ<double>;
template class BlockSparseMatrixEZ<float>;

DEAL_II_NAMESPACE_CLOSE
