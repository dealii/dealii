// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2018 by the deal.II authors
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

// test calling MultithreadInfo::set_thread_limit() multiple times and the
// effect for running tasks

#include <deal.II/base/thread_management.h>

#include "../tests.h"

Threads::Mutex mutex;
unsigned int   n_running;
unsigned int   n_max_running;

void
a_task()
{
  mutex.lock();
  n_running++;
  n_max_running = std::max(n_max_running, n_running);
  mutex.unlock();

  // Sleep some time to make sure all other concurrently running tasks enter
  // here.
  std::this_thread::sleep_for(std::chrono::seconds(1));

  mutex.lock();
  n_running--;
  mutex.unlock();
}


void
test()
{
  Threads::TaskGroup<> tg;

  n_running     = 0;
  n_max_running = 0;

  // force all tasks to wait until we are done starting
  mutex.lock();

  for (unsigned int t = 0; t < 10; ++t)
    tg += Threads::new_task(a_task);

  // now let the tasks run
  mutex.unlock();

  tg.join_all();
}


int
main()
{
  initlog();

  const unsigned int n = testing_max_num_threads();

  // if we have a machine with less than 5 cores fake the output:
  if (MultithreadInfo::n_threads() < n)
    {
      // you should buy more cores, seriously.
      deallog << "* default thread limit for tests: " << n << std::endl;
      deallog << "max concurrent running: " << n << " should be: " << n
              << std::endl;
    }
  else if (MultithreadInfo::n_threads() == n)
    {
      deallog << "* default thread limit for tests: "
              << MultithreadInfo::n_threads() << std::endl;
      test();
    }
  else
    Assert(false,
           ExcInternalError(
             "did somebody mess with LimitConcurrency in tests.h?"));

  deallog << "* now with thread limit 1:" << std::endl;
  MultithreadInfo::set_thread_limit(1);
  test();
  deallog << "* now with thread limit 2:" << std::endl;
  MultithreadInfo::set_thread_limit(2);
  test();

  deallog << "OK" << std::endl;
}
