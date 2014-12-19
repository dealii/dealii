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



// BlockSparsityPattern::column_number is broken

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/block_sparsity_pattern.h>

#include <iomanip>
#include <fstream>

int main()
{
  initlog();

  std::vector<types::global_dof_index> row_blocks(2);
  row_blocks[0] = 10;
  row_blocks[1] = 5;

  BlockCompressedSimpleSparsityPattern csp (row_blocks, row_blocks);

  csp.reinit(2,2);
  csp.block(0,0).reinit(10,10);
  csp.block(0,1).reinit(10,5);
  csp.block(1,0).reinit(5,10);
  csp.block(1,1).reinit(5,5);
  csp.collect_sizes();
  csp.add(5,5);
  csp.add(12,12);
  csp.add(5,11);
  csp.add(11,3);

  csp.print(deallog.get_file_stream());


  for (types::global_dof_index row=0; row<csp.n_rows(); ++row)
    {
      types::global_dof_index rlen = csp.row_length(row);

      for (types::global_dof_index c=0; c<rlen; ++c)
        {
          types::global_dof_index column = csp.column_number(row, c);
          deallog << row << "," << column << std::endl;
        }
    }

  deallog << "OK" << std::endl;
}
