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
// here: ComponentMask::represents_the_all_selected_mask


#include "../tests.h"
#include <deal.II/fe/component_mask.h>

#include <fstream>
#include <iomanip>




void test ()
{
  // test for an initialized mask
  Assert (ComponentMask(12,false).represents_the_all_selected_mask() == false,
          ExcInternalError());
  // note the semantics of the following as
  // described in the documentation
  Assert (ComponentMask(12,true).represents_the_all_selected_mask() == false,
          ExcInternalError());
  // test for an empty mask
  Assert (ComponentMask().represents_the_all_selected_mask() == true,
          ExcInternalError());

  deallog << "OK" << std::endl;
}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test();
}
