// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// start several tasks and see that return_value() waits for the correct
// result
//
// this test can also be used nicely to check that scheduling tasks works
// correctly. for example, on my 2-core laptop, it produces a load of 2
// (i.e. 2 tasks running concurrently) and runs with CPU time twice the wall
// time. on down, an 8-core machine, the load is 8 and the CPU time is 8 times
// the wall time

#include <deal.II/base/thread_management.h>

#include "../tests.h"


int
test(int i)
{
  int k = 1;
  for (unsigned int j = 0; j < 100000000; ++j)
    k += j % 17 + i;
  return k;
}



int
main()
{
  initlog();

  Threads::Task<int> tasks[] = {Threads::new_task(test, 1),
                                Threads::new_task(test, 2),
                                Threads::new_task(test, 3),
                                Threads::new_task(test, 4),
                                Threads::new_task(test, 5),
                                Threads::new_task(test, 6),
                                Threads::new_task(test, 7),
                                Threads::new_task(test, 8)};

  for (unsigned int i = 0; i < sizeof(tasks) / sizeof(tasks[0]); ++i)
    deallog << i << ' ' << tasks[i].return_value() << std::endl;

  deallog << "OK" << std::endl;
}
