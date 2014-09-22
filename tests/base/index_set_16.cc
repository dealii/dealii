// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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


// test IndexSet::subtract_set further

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdlib.h>

#include <deal.II/base/index_set.h>


void test ()
{
  IndexSet is1 (100);
  IndexSet is2 (100);

  is1.add_range(0,10);
  is1.add_range(20,100);

  is2.add_range(0,50);
//  is2.add_range(10,15);

  IndexSet is3 = is1;
  is3.subtract_set(is2);

  is1.print(deallog);
  is2.print(deallog);
  is3.print(deallog);

  for (unsigned int i=0; i<is3.size(); ++i)
    {
      Assert ((is1.is_element(i) && !is2.is_element(i))
              ==
              is3.is_element(i),
              ExcInternalError());
    }

  deallog << is3.index_within_set(51) << std::endl;
  Assert(is3.index_within_set(51)==1, ExcInternalError());

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
