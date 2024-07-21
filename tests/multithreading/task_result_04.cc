// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test TaskResult::try_emplace_task(). Make sure that any second call
// to this function does not actually change anything.

#include <deal.II/base/task_result.h>

#include "../tests.h"



int
main()
{
  initlog();

  // Create a default-constructed TaskResult object not associated with a task:
  Threads::TaskResult<int> task_result;

  // Next run 10 threads in parallel that all try to set a task
  // object. Calling TaskResult::operator=() will not work because
  // some thread will try to set a task while there is already a
  // running task, and that is considered an error.
  //
  // On the other hand, exactly one of the threads will succeed in
  // emplacing a task with try_emplace_task().
  //
  // To make things a bit more interesting, and make sure that there
  // really is contention, we let all the threads spin until a
  // condition variable changes state. In other words, they will all
  // try *at more or less the same time* to call try_emplace_task().
  std::atomic<bool> condition(false);

  const unsigned int     n_threads = 10;
  std::list<std::thread> threads;
  for (unsigned int t = 0; t < n_threads; ++t)
    threads.emplace_back(std::thread([&]() {
      while (condition.load() == false)
        /* keep checking */;

      task_result.try_emplace_task([]() { return 42; });

      // Now also try to emplace a task that will have a different
      // result: 12 instead of 42. None of these try_emplace_back()
      // calls will succeed because at the time we reach them, at
      // least one of the previous calls to the function returning 42
      // will have succeeded (on this or another thread). As a
      // consequence, the task result will have to be 42, regardless
      // of what we set here.
      task_result.try_emplace_task([]() { return 12; });
    }));

  // Now let all of these threads go to work:
  condition = true;

  // Wait for all of these threads to finish trying to emplace a task:
  for (auto &t : threads)
    t.join();

  // Wait for that background task to finish, then get its
  // value:
  deallog << task_result.value() << std::endl;

  // Ensure that we can continue to ask for the return value:
  deallog << task_result.value() << std::endl;
}
