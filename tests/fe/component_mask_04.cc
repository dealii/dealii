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


// tests for the ComponentMask class
//
// here: test ComponentMask::size()


#include "../tests.h"
#include <deal.II/fe/component_mask.h>

#include <fstream>
#include <iomanip>




void test ()
{
  Assert (ComponentMask(12, false).size() == 12, ExcInternalError());
  Assert (ComponentMask().size() == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}


int main()
{
  std::ofstream logfile ("component_mask_04/output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}
