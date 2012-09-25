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
// here: test that creating a mask that's empty is always true


#include "../tests.h"
#include <deal.II/fe/block_mask.h>

#include <fstream>
#include <iomanip>




void test ()
{
  BlockMask m;
  Assert (m[0] == true, ExcInternalError());
  Assert (m[42] == true, ExcInternalError());
  Assert (m[1000000000] == true, ExcInternalError());

  deallog << "OK" << std::endl;
}


int main()
{
  std::ofstream logfile ("block_mask_01/output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}
