// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2022 by the deal.II authors
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


// Check TaskGroup::return_values()


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

  // Create a bunch of tasks, and put them into a TaskGroup
  Threads::TaskGroup<int> tg;
  for (unsigned int i = 0; i < 10; ++i)
    tg += Threads::new_task(&test, i);

  // Then wait for them all to finish and output their respective
  // return values. They should be ordered the same way the tasks were
  // created, and should be twice the input argument:
  for (const int r : tg.return_values())
    deallog << r << std::endl;
}
