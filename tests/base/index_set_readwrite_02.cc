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


// test IndexSet::block_read() and block_write()

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <stdlib.h>

#include <deal.II/base/index_set.h>


void test ()
{
  IndexSet is1 (100);

  is1.add_range(0,10);
  is1.add_range(20,100);

  {
    std::ofstream out("a.idxset");
    is1.block_write(out);
  }

  IndexSet is2;
  std::ifstream in("a.idxset");
  is2.block_read(in);

  Assert(is1 == is2, ExcInternalError());


  IndexSet is3(11);
  is3.add_range(3,5);
  std::ifstream in2("a.idxset");
  is3.block_read(in2);

  Assert(is1 == is3, ExcInternalError());

  deallog << "OK" << std::endl;

  std::remove ("a.idxset");
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
