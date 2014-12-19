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



// tests for the BlockMask class
//
// here: BlockMask::n_selected_blocks


#include "../tests.h"
#include <deal.II/fe/block_mask.h>

#include <fstream>
#include <iomanip>




void test ()
{
  deal_II_exceptions::disable_abort_on_exception();

  // test for an initialized mask
  Assert (BlockMask(12,true).n_selected_blocks() == 12,
          ExcInternalError());
  Assert (BlockMask(12,true).n_selected_blocks(12) == 12,
          ExcInternalError());
  // test for an empty mask
  Assert (BlockMask().n_selected_blocks(12) == 12,
          ExcInternalError());
  Assert (BlockMask().n_selected_blocks(13) == 13,
          ExcInternalError());


  deallog << "OK" << std::endl;

  // this now must throw an exception,
  // though:
  try
    {
      Assert (BlockMask(12,true).n_selected_blocks(13) == 12,
              ExcInternalError());
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }

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
