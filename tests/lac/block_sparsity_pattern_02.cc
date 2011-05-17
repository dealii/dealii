//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


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
  std::ofstream logfile("block_sparsity_pattern_02/output");
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  std::vector<unsigned int> row_blocks(4);
  row_blocks[0] = 4;
  row_blocks[1] = 5;
  row_blocks[2] = 1;
  row_blocks[3] = 4;
  std::vector<unsigned int> col_blocks(3);
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
