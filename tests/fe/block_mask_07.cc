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
// here: BlockMask::first_selected_block


#include <deal.II/fe/block_mask.h>

#include "../tests.h"



void
test()
{
  std::vector<bool> v(12, false);
  v[3] = true;

  BlockMask m(v);

  // test for an initialized mask
  Assert(m.first_selected_block() == 3, ExcInternalError());
  Assert(BlockMask(12, true).first_selected_block() == 0, ExcInternalError());
  // test for an empty mask
  Assert(BlockMask().first_selected_block(12) == 0, ExcInternalError());

  deallog << "OK" << std::endl;

  // the following should yield an exception:
  try
    {
      Assert(BlockMask(12, true).first_selected_block(13) == 0,
             ExcInternalError());
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }

  // as should this:
  try
    {
      Assert(BlockMask(12, false).first_selected_block() == 0,
             ExcInternalError());
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
}


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();
  deallog << std::setprecision(4);

  test();
}
