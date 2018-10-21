// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// tests for the ComponentMask class
//
// here: test that creating a mask that's empty is always true


#include <deal.II/fe/component_mask.h>

#include "../tests.h"



void
test()
{
  ComponentMask m;
  DEAL_II_AssertThrow(m[0] == true, ExcInternalError());
  DEAL_II_AssertThrow(m[42] == true, ExcInternalError());
  DEAL_II_AssertThrow(m[1000000000] == true, ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(4);

  deallog.attach(logfile);

  test();
}
