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
// here: ComponentMask::operator==


#include "../tests.h"
#include <deal.II/fe/component_mask.h>

#include <fstream>
#include <iomanip>




void test ()
{
  std::vector<bool> v1(12);
  for (unsigned int i=0; i<v1.size(); ++i)
    v1[i] = (i%3 == 0);
  std::vector<bool> v2(12);
  for (unsigned int i=0; i<v2.size(); ++i)
    v2[i] = (i%4 == 0);

  std::vector<bool> v(12);
  for (unsigned int i=0; i<v.size(); ++i)
    v[i] = (v1[i] || v2[i]);

  ComponentMask m1(v1);
  ComponentMask m2(v2);
  ComponentMask m = m1 | m2;

  // verify equality
  Assert (m == ComponentMask(v),
          ExcInternalError());
  Assert (!(m == m1),
          ExcInternalError());
  Assert (!(m == ComponentMask(v1)),
          ExcInternalError());
  Assert (!(m == ComponentMask(v2)),
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
