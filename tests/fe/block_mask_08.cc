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
// here: BlockMask::represents_the_all_selected_mask


#include "../tests.h"
#include <deal.II/fe/block_mask.h>

#include <fstream>
#include <iomanip>




void test ()
{
				   // test for an initialized mask
  Assert (BlockMask(12,false).represents_the_all_selected_mask() == false,
	  ExcInternalError());
				   // note the semantics of the following as
				   // described in the documentation
  Assert (BlockMask(12,true).represents_the_all_selected_mask() == false,
	  ExcInternalError());
				   // test for an empty mask
  Assert (BlockMask().represents_the_all_selected_mask() == true,
	  ExcInternalError());

  deallog << "OK" << std::endl;
}


int main()
{
  std::ofstream logfile ("block_mask_08/output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}
