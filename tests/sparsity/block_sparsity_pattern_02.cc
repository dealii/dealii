// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// BlockSparsityPattern::copy_from only copied n_block_rows x
// n_block_rows blocks, forgetting any additional columns there may
// have been (or trying to copy columns that aren't there


#include <deal.II/lac/block_sparsity_pattern.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(2) << std::fixed;

  std::vector<types::global_dof_index> row_blocks(4);
  row_blocks[0] = 4;
  row_blocks[1] = 5;
  row_blocks[2] = 1;
  row_blocks[3] = 4;
  std::vector<types::global_dof_index> col_blocks(3);
  col_blocks[0] = 2;
  col_blocks[1] = 3;
  col_blocks[2] = 2;

  BlockDynamicSparsityPattern bcsp(row_blocks, col_blocks);

  BlockSparsityPattern bsp;

  // we used to access an invalid
  // block number here
  bsp.copy_from(bcsp);

  deallog << "OK" << std::endl;
}
