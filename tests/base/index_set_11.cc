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


// test IndexSet::add_indices

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/index_set.h>


void test ()
{
  IndexSet index_set (10);

  {
    const unsigned int array[] = { 2, 3, 4 };
    index_set.add_indices ((const unsigned int *)array, array+sizeof(array)/sizeof(array[0]));
  }
  {
    const unsigned int array[] = { 7,6 };
    index_set.add_indices ((const unsigned int *)array, array+sizeof(array)/sizeof(array[0]));
  }

  Assert (index_set.is_contiguous() == false, ExcInternalError());

  for (unsigned int i=0; i<index_set.size(); ++i)
    deallog << i << ' ' << (index_set.is_element(i) ? "true" : "false")
            << std::endl;


  {
    const unsigned int array[] = { 5 };
    index_set.add_indices ((const unsigned int *)array, array+sizeof(array)/sizeof(array[0]));
  }

  Assert (index_set.is_contiguous() == true, ExcInternalError());

  for (unsigned int i=0; i<index_set.size(); ++i)
    deallog << i << ' ' << (index_set.is_element(i) ? "true" : "false")
            << std::endl;


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
