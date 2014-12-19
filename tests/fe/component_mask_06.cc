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
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }

}


int main()
{
  deal_II_exceptions::disable_abort_on_exception();

  std::ofstream logfile ("output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}
