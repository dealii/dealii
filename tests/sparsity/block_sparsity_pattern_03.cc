// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// BlockSparsityPattern::column_number is broken

#include <deal.II/lac/block_sparsity_pattern.h>

#include "../tests.h"


int
main()
{
  initlog();

  std::vector<types::global_dof_index> row_blocks(2);
  row_blocks[0] = 10;
  row_blocks[1] = 5;

  BlockDynamicSparsityPattern csp(row_blocks, row_blocks);

  csp.reinit(2, 2);
  csp.block(0, 0).reinit(10, 10);
  csp.block(0, 1).reinit(10, 5);
  csp.block(1, 0).reinit(5, 10);
  csp.block(1, 1).reinit(5, 5);
  csp.collect_sizes();
  csp.add(5, 5);
  csp.add(12, 12);
  csp.add(5, 11);
  csp.add(11, 3);

  csp.print(deallog.get_file_stream());


  for (types::global_dof_index row = 0; row < csp.n_rows(); ++row)
    {
      types::global_dof_index rlen = csp.row_length(row);

      for (types::global_dof_index c = 0; c < rlen; ++c)
        {
          types::global_dof_index column = csp.column_number(row, c);
          deallog << row << ',' << column << std::endl;
        }
    }

  deallog << "OK" << std::endl;
}
