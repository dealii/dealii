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


// start several tasks and see that return_value() waits for the correct
// result
//
// this test can also be used nicely to check that scheduling tasks works
// correctly. for example, on my 2-core laptop, it produces a load of 2
// (i.e. 2 tasks running concurrently) and runs with CPU time twice the wall
// time. on down, an 8-core machine, the load is 8 and the CPU time is 8 times
// the wall time

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <unistd.h>

#include <deal.II/base/thread_management.h>


int test (int i)
{
  int k = 1;
  for (unsigned int j=0; j<100000000; ++j)
    k += j%17+i;
  return k;
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Threads::Task<int> tasks[] =
  {
    Threads::new_task (test, 1),
    Threads::new_task (test, 2),
    Threads::new_task (test, 3),
    Threads::new_task (test, 4),
    Threads::new_task (test, 5),
    Threads::new_task (test, 6),
    Threads::new_task (test, 7),
    Threads::new_task (test, 8)
  };

  for (unsigned int i=0; i<sizeof(tasks)/sizeof(tasks[0]); ++i)
    deallog << i << ' ' << tasks[i].return_value() << std::endl;

  deallog << "OK" << std::endl;
}
