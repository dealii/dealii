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
// test checks the move constructor.

#include <deal.II/base/task_result.h>

#include "../tests.h"



int
main()
{
  initlog();

  // Create a background task that sleeps for a while and then
  // issues its result. Associate this background task with
  // a TaskResult object.
  Threads::TaskResult<int> t = Threads::new_task([]() {
    std::this_thread::sleep_for(std::chrono::seconds(1));
    return 42;
  });

  // Move the task
  Threads::TaskResult<int> tt = std::move(t);

  // Wait for that background task to finish, then get its
  // value:
  deallog << tt.value() << std::endl;

  // Ensure that we can continue to ask for the return value:
  deallog << tt.value() << std::endl;
}
