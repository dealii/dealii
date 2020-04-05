// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2018 by the deal.II authors
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



// tests for the BlockMask class
//
// here: BlockMask::represents_n_blocks


#include <deal.II/fe/block_mask.h>

#include "../tests.h"



void
test()
{
  // test for an initialized mask
  AssertThrow(BlockMask(12, false).represents_n_blocks(12) == true,
              ExcInternalError());
  AssertThrow(BlockMask(12, false).represents_n_blocks(13) == false,
              ExcInternalError());
  // test for an empty mask
  AssertThrow(BlockMask().represents_n_blocks(12) == true, ExcInternalError());
  AssertThrow(BlockMask().represents_n_blocks(13) == true, ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test();
}
