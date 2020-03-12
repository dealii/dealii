// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2018 by the deal.II authors
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

  deallog << "blocks: " << csp.n_block_rows() << "x" << csp.n_block_cols()
          << std::endl;
  deallog << "size: " << csp.n_rows() << "x" << csp.n_cols() << std::endl;
  deallog << "size block(1,0):" << csp.block(1, 0).n_rows() << "x"
          << csp.block(1, 0).n_cols() << std::endl;

  csp.print(deallog.get_file_stream());

  deallog << std::endl;
  for (unsigned int i = 0; i < 3; ++i)
    {
      deallog << "row " << i << ": ";

      for (unsigned int j = 0; j < csp.row_length(i); ++j)
        deallog << csp.column_number(i, j) << " ";

      deallog << std::endl;
    }

  return 0;
}
