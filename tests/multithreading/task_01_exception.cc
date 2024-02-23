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


// Like the _01 test, but check that if the function called on the
// task throws an exception, that Task::wait() continues to work but
// Task::return_value() throws the exception.

#include <deal.II/base/thread_management.h>

#include "../tests.h"


int
test()
{
  std::this_thread::sleep_for(std::chrono::seconds(3));

  throw 42;
}


int
main()
{
  initlog();

  Threads::Task<int> t = Threads::new_task(test);

  // Here, an exception should be triggered:
  try
    {
      t.join();
    }
  catch (int i)
    {
      Assert(i == 42, ExcInternalError());
      deallog << "OK" << std::endl;
    }
}
