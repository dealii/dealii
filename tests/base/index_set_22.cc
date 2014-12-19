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


// test IndexSet::fill_binary_vector

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

  // randomly add 90 elements to each
  // set, some of which may be
  // repetitions of previous ones
  for (unsigned int i=0; i<9*is1.size()/10; ++i)
    is1.add_index (Testing::rand() % is1.size());

  std::vector<bool> zeros_and_ones(is1.size());
  is1.fill_binary_vector (zeros_and_ones);

  deallog << "Original index set: " << std::endl;
  is1.print(deallog);

  for (unsigned int i=0; i<is1.size(); i++)
    Assert(is1.is_element(i) == zeros_and_ones[i],
           ExcInternalError());

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
