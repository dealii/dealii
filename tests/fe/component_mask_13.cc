// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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



// tests for the ComponentMask class
//
// here: test conversion from component mask to block mask and back. start with
// blocks where the _12 test starts with components, and use an element that
// does have fewer blocks than components


#include "../tests.h"
#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_raviart_thomas.h>

#include <fstream>
#include <iomanip>




void test ()
{
  // the following element has more
  // components than blocks
  FESystem<2> fe (FE_RaviartThomas<2>(1), 3);

  // try all possible block
  // masks, which we encode as bit
  // strings
  for (unsigned int int_mask=0; int_mask<(1U<<fe.n_blocks()); ++int_mask)
    {
      std::vector<bool> block_mask (fe.n_blocks());
      for (unsigned int c=0; c<fe.n_blocks(); ++c)
        block_mask[c] = (int_mask & (1<<c));

      // make sure that the round-trip works
      Assert (BlockMask(block_mask) == fe.block_mask(fe.component_mask(BlockMask(block_mask))),
              ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}
