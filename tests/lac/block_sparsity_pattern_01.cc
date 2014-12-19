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



// Test reinit(BlockIndices...)


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

  BlockSparsityPattern sparsity;
  std::vector<types::global_dof_index> row_blocks(4);
  row_blocks[0] = 4;
  row_blocks[1] = 5;
  row_blocks[2] = 1;
  row_blocks[3] = 4;
  std::vector<types::global_dof_index> col_blocks(3);
  col_blocks[0] = 2;
  col_blocks[1] = 3;
  col_blocks[2] = 2;
  BlockIndices rows(row_blocks);
  BlockIndices cols(col_blocks);

  std::vector<std::vector<unsigned int> >
  row_length(cols.size(),
             std::vector<unsigned int>(rows.total_size()));
  for (unsigned int jb=0; jb<row_length.size(); ++jb)
    for (unsigned int i=0; i<row_length[jb].size(); ++i)
      {
        const unsigned int d = col_blocks[jb]-1;
        row_length[jb][i] = (i+1) % d +1;
      }

  for (unsigned int j=0; j<row_length.size(); ++j)
    {
      for (unsigned int i=0; i<row_length[j].size(); ++i)
        deallog << ' ' << row_length[j][i];
      deallog << std::endl;
    }

  sparsity.reinit(rows, cols, row_length);

  for (unsigned int ib = 0; ib<rows.size(); ++ib)
    for (unsigned int i=0; i<rows.block_size(ib); ++i)
      {
        const unsigned int ii = rows.local_to_global(ib,i);
        for (unsigned int jb = 0; jb<cols.size(); ++jb)
          for (unsigned int j=0; j<row_length[jb][ii]; ++j)
            sparsity.add(ii, cols.local_to_global(jb,j));
      }

  sparsity.print(logfile);

  return 0;
}
