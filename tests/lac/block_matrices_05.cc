// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test that the frobenius_norm function returns the correct value
// for a block matrix.


#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>

#include <algorithm>

#include "../tests.h"



int
main()
{
  initlog();

  BlockSparsityPattern bsp(2, 2);
  // set sizes
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      bsp.block(i, j).reinit(5, 5, 5);
  bsp.collect_sizes();

  // make a full matrix
  for (unsigned int row = 0; row < 10; ++row)
    for (unsigned int i = 0; i < 10; ++i)
      bsp.add(row, i);
  bsp.compress();

  BlockSparseMatrix<double> bsm(bsp);

  for (unsigned int i = 0; i < bsm.m(); ++i)
    for (unsigned int j = 0; j < bsm.n(); ++j)
      {
        bsm.add(i, j, 1.0);
      }

  const double accuracy          = 1e-12;
  const double correct_frob_norm = 10.0;
  const double frob_norm         = bsm.frobenius_norm();

  Assert(std::fabs(frob_norm - correct_frob_norm) < accuracy,
         ExcInternalError());

  deallog << "OK" << std::endl;

  return 0;
}
