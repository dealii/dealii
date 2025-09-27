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


// A test for basic functionality of the TaskResult function. This
// test checks TaskResult::empty()'s documented functionality.

#include <deal.II/base/task_result.h>

#include "../tests.h"



int
main()
{
  initlog();

  {
    Threads::TaskResult<int> t;
    AssertThrow(t.empty() == true, ExcInternalError());
  }


  {
    // Create a background task that sleeps for a while and then
    // issues its result. Associate this background task with
    // a TaskResult object.
    Threads::TaskResult<int> t = Threads::new_task([]() {
      std::this_thread::sleep_for(std::chrono::seconds(1));
      return 42;
    });
    AssertThrow(t.empty() == false, ExcInternalError());

    t.join();
    AssertThrow(t.empty() == false, ExcInternalError());
  }

  // Same idea, but move construct the task
  {
    // Create a background task that sleeps for a while and then
    // issues its result. Associate this background task with
    // a TaskResult object.
    Threads::TaskResult<int> t = Threads::new_task([]() {
      std::this_thread::sleep_for(std::chrono::seconds(1));
      return 42;
    });

    Threads::TaskResult<int> tt = std::move(t);
    AssertThrow(t.empty() == true, ExcInternalError());
    AssertThrow(tt.empty() == false, ExcInternalError());

    // Ensure the task has finished before destroying the TaskResult
    // object:
    tt.join();
  }

  {
    // TaskResult::clear() results in an empty object, but we have to
    // wait in order to call clear(). We can do this by explicitly
    // asking for the resulting value...
    Threads::TaskResult<int> t     = Threads::new_task([]() {
      std::this_thread::sleep_for(std::chrono::seconds(1));
      return 42;
    });
    const int                value = t.value();
    AssertThrow(value == 42, ExcInternalError());
    AssertThrow(t.empty() == false, ExcInternalError());

    t.clear();
    AssertThrow(t.empty() == true, ExcInternalError());
  }


  {
    // ...or using the weaker t.join() call.
    Threads::TaskResult<int> t = Threads::new_task([]() {
      std::this_thread::sleep_for(std::chrono::seconds(1));
      return 42;
    });
    t.join();
    AssertThrow(t.empty() == false, ExcInternalError());

    t.clear();
    AssertThrow(t.empty() == true, ExcInternalError());
  }

  deallog << "OK" << std::endl;
}
