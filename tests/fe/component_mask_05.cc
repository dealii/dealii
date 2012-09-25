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
// here: ComponentMask::represents_n_components


#include "../tests.h"
#include <deal.II/fe/component_mask.h>

#include <fstream>
#include <iomanip>




void test ()
{
				   // test for an initialized mask
  Assert (ComponentMask(12,false).represents_n_components(12) == true,
	  ExcInternalError());
  Assert (ComponentMask(12,false).represents_n_components(13) == false,
	  ExcInternalError());
				   // test for an empty mask
  Assert (ComponentMask().represents_n_components(12) == true,
	  ExcInternalError());
  Assert (ComponentMask().represents_n_components(13) == true,
	  ExcInternalError());

  deallog << "OK" << std::endl;
}


int main()
{
  std::ofstream logfile ("component_mask_05/output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}
