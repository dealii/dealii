//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


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
  std::ofstream logfile ("block_mask_06/output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}
