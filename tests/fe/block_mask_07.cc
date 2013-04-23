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
// here: BlockMask::first_selected_block


#include "../tests.h"
#include <deal.II/fe/block_mask.h>

#include <fstream>
#include <iomanip>




void test ()
{
  std::vector<bool> v(12,false);
  v[3] = true;

  BlockMask m(v);

				   // test for an initialized mask
  Assert (m.first_selected_block() == 3,
	  ExcInternalError());
  Assert (BlockMask(12,true).first_selected_block() == 0,
	  ExcInternalError());
				   // test for an empty mask
  Assert (BlockMask().first_selected_block(12) == 0,
	  ExcInternalError());

  deallog << "OK" << std::endl;

				   // the following should yield an exception:
  try 
    {
      Assert (BlockMask(12,true).first_selected_block(13) == 0,
	      ExcInternalError());
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
  
				   // as should this: 
  try 
    {
      Assert (BlockMask(12,false).first_selected_block() == 0,
	      ExcInternalError());
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
}


int main()
{
  deal_II_exceptions::disable_abort_on_exception();

  std::ofstream logfile ("block_mask_07/output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}
