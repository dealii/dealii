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


// a simple test for task based programming: just create a task and wait for
// its completion. the task sleeps for a bit to make sure that the waiting
// code works alright

#include <deal.II/base/thread_management.h>

#include "../tests.h"


void
test()
{
  std::this_thread::sleep_for(std::chrono::seconds(3));
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  Threads::Task<> t = Threads::new_task(test);
  t.join();
}
