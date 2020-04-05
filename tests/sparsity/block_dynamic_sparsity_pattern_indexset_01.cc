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



// Test BlockDynamicSparsityPattern with IndexSets


#include <deal.II/lac/block_sparsity_pattern.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(2) << std::fixed;

  IndexSet a(5);
  IndexSet b(3);
  a.add_index(0);
  a.add_index(3);
  b.add_index(0);
  std::vector<IndexSet> part;
  part.push_back(a);
  part.push_back(b);
  part.push_back(a);
  BlockDynamicSparsityPattern csp(part);

  deallog << "blocks: " << csp.n_block_rows() << "x" << csp.n_block_cols()
          << std::endl;
  deallog << "size: " << csp.n_rows() << "x" << csp.n_cols() << std::endl;
  deallog << "size block(1,0):" << csp.block(1, 0).n_rows() << "x"
          << csp.block(1, 0).n_cols() << std::endl;

  part.pop_back();
  csp.reinit(part); // also check the reinit variant.

  deallog << "blocks: " << csp.n_block_rows() << "x" << csp.n_block_cols()
          << std::endl;
  deallog << "size: " << csp.n_rows() << "x" << csp.n_cols() << std::endl;
  deallog << "size block(1,0):" << csp.block(1, 0).n_rows() << "x"
          << csp.block(1, 0).n_cols() << std::endl;

  for (int i = 0; i < 13; ++i)
    {
      if (i == 0 || i == 3 || i == 5)
        {
          csp.add(i, 0);
          csp.add(i, i);
        }
    }

  csp.print(deallog.get_file_stream());

  return 0;
}
