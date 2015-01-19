// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/lac/block_sparse_matrix_ez.h>
#include <deal.II/lac/block_sparse_matrix_ez.templates.h>

DEAL_II_NAMESPACE_OPEN

// explicit instantiations
template class BlockSparseMatrixEZ<double>;
template class BlockSparseMatrixEZ<float>;

DEAL_II_NAMESPACE_CLOSE
