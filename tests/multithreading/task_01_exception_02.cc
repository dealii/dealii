// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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
      t.return_value();
    }
  catch (int i)
    {
      Assert(i == 42, ExcInternalError());
      deallog << "OK" << std::endl;
    }
}
