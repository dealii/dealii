// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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


// test IndexSet::operator &

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/index_set.h>


void test ()
{
  IndexSet is1 (20);
  IndexSet is2 (20);

  is1.add_index (2);
  is1.add_index (3);
  is1.add_index (4);
  is1.add_index (6);
  is1.add_index (7);

  is2.add_range (4,9);

  IndexSet is3 = is1 & is2;

  for (unsigned int i=0; i<is3.size(); ++i)
    {
      deallog << i << ' ' << (is3.is_element(i) ? "true" : "false")
              << std::endl;

      Assert ((is1.is_element(i) && is2.is_element(i))
              ==
              is3.is_element(i),
              ExcInternalError());
    }

  // some sanity tests
  Assert ((is1 & is2) == (is2 & is1), ExcInternalError());
  Assert ((is1 & is3) == (is2 & is3), ExcInternalError());
  Assert ((is1 & is3) == is3, ExcInternalError());
  Assert ((is3 & is1) == is3, ExcInternalError());
  Assert ((is3 & is3) == is3, ExcInternalError());

  deallog << "OK" << std::endl;
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
