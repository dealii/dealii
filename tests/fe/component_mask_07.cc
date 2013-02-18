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
// here: ComponentMask::first_selected_component


#include "../tests.h"
#include <deal.II/fe/component_mask.h>

#include <fstream>
#include <iomanip>




void test ()
{
  std::vector<bool> v(12,false);
  v[3] = true;

  ComponentMask m(v);

				   // test for an initialized mask
  Assert (m.first_selected_component() == 3,
	  ExcInternalError());
  Assert (ComponentMask(12,true).first_selected_component() == 0,
	  ExcInternalError());
				   // test for an empty mask
  Assert (ComponentMask().first_selected_component(12) == 0,
	  ExcInternalError());

  deallog << "OK" << std::endl;

				   // the following should yield an exception:
  try 
    {
      Assert (ComponentMask(12,true).first_selected_component(13) == 0,
	      ExcInternalError());
    }
  catch (...)
    {
    }
  
				   // as should this:
  try
    {
      Assert (ComponentMask(12,false).first_selected_component() == 0,
	      ExcInternalError());
    }
  catch (...)
    {
    }
}


int main()
{
  deal_II_exceptions::disable_abort_on_exception();

  std::ofstream logfile ("component_mask_07/output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}
