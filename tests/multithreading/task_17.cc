// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Allow only 2 threads (= concurrently running tasks) but create 3
// and make sure each can wait for the next to finish. This requires
// that in the underlying implementation, when we wait for a task we
// go back into the scheduler and execute more tasks. In other words,
// waiting for another task doesn't just block the current task (which
// would lead to deadlocks because we're already running the maximum
// number of current tasks) but gives the scheduler the opportunity to
// work on other tasks.


#include <deal.II/base/thread_management.h>

#include "../tests.h"

void
bottom()
{
  deallog << "      Starting task at the bottom" << std::endl;
  deallog << "        ... ... ..." << std::endl;
  std::this_thread::sleep_for(std::chrono::seconds(1));
  deallog << "      Ending task at the bottom" << std::endl;
}

void
middle()
{
  deallog << "    Starting task in the middle" << std::endl;
  auto t = Threads::new_task([]() { bottom(); });
  t.join();
  deallog << "    Ending task in the middle" << std::endl;
}

void
top()
{
  deallog << "  Starting task at the top" << std::endl;
  auto t = Threads::new_task([]() { middle(); });
  t.join();
  deallog << "  Ending task at the top" << std::endl;
}


int
main()
{
  initlog();

  MultithreadInfo::set_thread_limit(2);

  deallog << "Starting task in main()" << std::endl;
  auto t = Threads::new_task([]() { top(); });
  t.join();
  deallog << "Done in main" << std::endl;
}
