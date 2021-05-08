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


// like task_01, but with return value

#include <deal.II/base/thread_management.h>

#include "../tests.h"


int
test()
{
  std::this_thread::sleep_for(std::chrono::seconds(3));
  return 42;
}



int
main()
{
  initlog();

  Threads::Task<int> t = Threads::new_task(test);
  AssertThrow(t.return_value() == 42, ExcInternalError());

  deallog << "OK" << std::endl;
}
