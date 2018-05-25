// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2017 by the deal.II authors
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


// verify Threads::Task::joinable()

#include <deal.II/base/thread_management.h>

#include <unistd.h>

#include "../tests.h"


void
test()
{
  sleep(3);
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  Threads::Task<> t;
  AssertThrow(t.joinable() == false, ExcInternalError());

  t = Threads::new_task(test);
  AssertThrow(t.joinable() == true, ExcInternalError());
  t.join();
}
