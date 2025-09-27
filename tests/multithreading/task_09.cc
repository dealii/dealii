// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// we used to synchronize child tasks with the one that spawned it by
// using a mutex, but this could lead to deadlocks if there are more
// tasks than processors available, see the emails on the mailing
// lists in late nov 2010.
//
// this test verifies that we no longer deadlock

#include <deal.II/base/multithread_info.h>
#include <deal.II/base/thread_management.h>

#include <chrono>
#include <thread>

#include "../tests.h"


void
test()
{
  std::this_thread::sleep_for(std::chrono::seconds(1));
}


void
outer()
{
  // wait for some time to make sure that really all outer tasks have been
  // started, then start a bunch of new tasks and wait for them to finish. it
  // used to be that the join() function used a mutex which could only be
  // acquired once the spawned task has finished, but since we already have so
  // many tasks running at this point, the children never get to run and so
  // never release the mutex that we try to acquire in join()
  std::this_thread::sleep_for(std::chrono::seconds(1));

  Threads::new_task(test).join();
}



int
main()
{
  initlog();

  {
    const std::size_t n_tasks = MultithreadInfo::n_cores() * 2;
#ifdef _POSIX_C_SOURCE
    // get rid of the environment variable if it exists (this may not be
    // possible under MSVC, but the test will still pass; it just won't test
    // what we want it to test)
    unsetenv("DEAL_II_NUM_THREADS");
#endif
    // Since this test relies on providing more tasks than there are cores,
    // reset the concurrency limit that exists in tests.h:
    MultithreadInfo::set_thread_limit(n_tasks);
    std::vector<Threads::Task<void>> tasks(n_tasks);
    for (std::size_t task_n = 0; task_n < n_tasks; ++task_n)
      {
        // We must save each task, otherwise it will join before the next loop
        // iteration
        tasks.emplace_back(Threads::new_task(outer));
      }
  }

  deallog << "OK" << std::endl;
}
