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
// here: ComponentMask::n_selected_components


#include "../tests.h"
#include <deal.II/fe/component_mask.h>

#include <fstream>
#include <iomanip>




void test ()
{
				   // test for an initialized mask
  Assert (ComponentMask(12,true).n_selected_components() == 12,
	  ExcInternalError());
  Assert (ComponentMask(12,true).n_selected_components(12) == 12,
	  ExcInternalError());
				   // test for an empty mask
  Assert (ComponentMask().n_selected_components(12) == 12,
	  ExcInternalError());
  Assert (ComponentMask().n_selected_components(13) == 13,
	  ExcInternalError());


  deallog << "OK" << std::endl;

				   // this now must throw an exception,
				   // though:
  try 
    {
      Assert (ComponentMask(12,true).n_selected_components(13) == 12,
	      ExcInternalError());
    }
  catch (...)
    {
    }
  
}


int main()
{
  deal_II_exceptions::disable_abort_on_exception();

  std::ofstream logfile ("component_mask_06/output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}
