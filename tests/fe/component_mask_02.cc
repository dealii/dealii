// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2014 by the deal.II authors
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
// here: test that creating a mask from a vector<bool> works


#include "../tests.h"
#include <deal.II/fe/component_mask.h>

#include <fstream>
#include <iomanip>




void test ()
{
  std::vector<bool> v(12);
  for (unsigned int i=0; i<v.size(); ++i)
    v[i] = (i%3 == 0);

  ComponentMask m(v);

  // verify equality
  for (unsigned int i=0; i<v.size(); ++i)
    Assert (m[i] == v[i], ExcInternalError());

  // this needs to throw an exception
  try
    {
      m[v.size()];
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }

  deallog << "OK" << std::endl;
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
