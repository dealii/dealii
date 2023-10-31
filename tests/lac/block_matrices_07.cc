// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2021 by the deal.II authors
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

// It was not possible to compare non-const iterators into block
// sparse matrices because we declared the wrong class as a
// friend. This test checks that this is fixed.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/block_sparse_matrix.h>

#include "../tests.h"

#include "../testmatrix.h"


dealii::BlockSparseMatrix<double> systemMatrix;
int
main()
{
  initlog();

  const unsigned int size = 5;
  unsigned int       dim  = (size - 1) * (size - 1);

  FDMatrix testproblem(size, size);

  BlockSparsityPattern block_structure(1, 1);
  block_structure.block(0, 0).reinit(dim, dim, 5);
  block_structure.collect_sizes();
  testproblem.five_point_structure(block_structure.block(0, 0));
  block_structure.compress();
  BlockSparseMatrix<double> block_A(block_structure);

  auto block_AIteratorEnd = block_A.end(0);

  auto block_AIteratorBegin = block_A.begin(0);

  if (block_AIteratorBegin != block_AIteratorEnd)
    deallog << "Not equal" << std::endl;

  deallog << "OK." << std::endl;
}
