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
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }

  // as should this:
  try
    {
      Assert (ComponentMask(12,false).first_selected_component() == 0,
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
