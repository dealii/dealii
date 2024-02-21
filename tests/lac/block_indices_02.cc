// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test BlockIndices::global_to_local() for very large number of blocks

#include <deal.II/lac/block_indices.h>

#include "../tests.h"



void
test()
{
  constexpr unsigned int            n_blocks   = 10000;
  constexpr types::global_dof_index block_size = 3;
  constexpr types::global_dof_index size       = n_blocks * block_size;
  BlockIndices                      idx(n_blocks, block_size);

  const std::vector<types::global_dof_index> is = {{0,
                                                    1,
                                                    block_size * 500 + 2,
                                                    block_size * 600,
                                                    block_size * 701 - 1,
                                                    size - 1}};

  deallog << "n_blocks:   " << n_blocks << std::endl
          << "block_size: " << block_size << std::endl
          << "size:       " << idx.size() << std::endl;
  for (auto i : is)
    {
      const auto p = idx.global_to_local(i);
      AssertDimension(p.first, i / block_size);
      AssertDimension(p.second, i % block_size);
      deallog << i << " -> (" << p.first << ',' << p.second << ')' << std::endl;
    }

  deallog << std::endl;
}


int
main()
{
  initlog();
  test();
}
