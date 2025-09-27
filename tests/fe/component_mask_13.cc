// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// tests for the ComponentMask class
//
// here: test conversion from component mask to block mask and back. start with
// blocks where the _12 test starts with components, and use an element that
// does have fewer blocks than components


#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include "../tests.h"



void
test()
{
  // the following element has more
  // components than blocks
  FESystem<2> fe(FE_RaviartThomas<2>(1), 3);

  // try all possible block
  // masks, which we encode as bit
  // strings
  for (unsigned int int_mask = 0; int_mask < (1U << fe.n_blocks()); ++int_mask)
    {
      std::vector<bool> block_mask(fe.n_blocks());
      for (unsigned int c = 0; c < fe.n_blocks(); ++c)
        block_mask[c] = (int_mask & (1 << c));

      // make sure that the round-trip works
      AssertThrow(BlockMask(block_mask) ==
                    fe.block_mask(fe.component_mask(BlockMask(block_mask))),
                  ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test();
}
