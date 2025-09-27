// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Like the _15 test, but check that arguments are properly copied.


#include <deal.II/base/thread_management.h>

#include "../tests.h"


// return a double, to make sure we correctly identify the return type
// of the expressions used in new_task(...)
int
test(int i)
{
  return 2 * i;
}



int
main()
{
  initlog();

  // Create a bunch of tasks, and put them into a TaskGroup. The test
  // also checks that the call to `new_task()` needs to copy the
  // passed in value for `i` rather than take a reference to it,
  // because it would otherwise use outdated values.
  Threads::TaskGroup<int> tg;
  for (unsigned int i = 0; i < 10; ++i)
    tg += Threads::new_task([](int ii) { return test(ii); }, i);

  // Then wait for them all to finish and output their respective
  // return values. They should be ordered the same way the tasks were
  // created, and should be twice the input argument:
  for (const int r : tg.return_values())
    deallog << r << std::endl;
}
