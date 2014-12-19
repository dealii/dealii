// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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



// BlockSparsityPattern::copy_from only copied n_block_rows x
// n_block_rows blocks, forgetting any additional columns there may
// have been (or trying to copy columns that aren't there


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/block_sparsity_pattern.h>

#include <iomanip>
#include <fstream>

int main()
{
  std::ofstream logfile("output");
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  std::vector<types::global_dof_index> row_blocks(4);
  row_blocks[0] = 4;
  row_blocks[1] = 5;
  row_blocks[2] = 1;
  row_blocks[3] = 4;
  std::vector<types::global_dof_index> col_blocks(3);
  col_blocks[0] = 2;
  col_blocks[1] = 3;
  col_blocks[2] = 2;

  BlockCompressedSparsityPattern bcsp (row_blocks, col_blocks);

  BlockSparsityPattern bsp;

  // we used to access an invalid
  // block number here
  bsp.copy_from (bcsp);

  deallog << "OK" << std::endl;
}
