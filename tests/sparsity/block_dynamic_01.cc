// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test BlockDynamicSparsityPattern column_number


#include <deal.II/lac/block_sparsity_pattern.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(2) << std::fixed;

  BlockDynamicSparsityPattern csp(2, 2);
  csp.block(0, 0).reinit(2, 2);
  csp.block(1, 0).reinit(1, 2);
  csp.block(0, 1).reinit(2, 1);
  csp.block(1, 1).reinit(1, 1);
  csp.collect_sizes();
  for (int i = 0; i < 3; ++i)
    csp.add(i, i);
  csp.add(1, 0);
  csp.add(1, 2);
  csp.compress();

  deallog << "blocks: " << csp.n_block_rows() << 'x' << csp.n_block_cols()
          << std::endl;
  deallog << "size: " << csp.n_rows() << 'x' << csp.n_cols() << std::endl;
  deallog << "size block(1,0):" << csp.block(1, 0).n_rows() << 'x'
          << csp.block(1, 0).n_cols() << std::endl;

  csp.print(deallog.get_file_stream());

  deallog << std::endl;
  for (unsigned int i = 0; i < 3; ++i)
    {
      deallog << "row " << i << ": ";

      for (unsigned int j = 0; j < csp.row_length(i); ++j)
        deallog << csp.column_number(i, j) << ' ';

      deallog << std::endl;
    }

  return 0;
}
