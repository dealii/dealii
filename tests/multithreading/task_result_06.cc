// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2009 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Test TaskResult::emplace_object().

#include <deal.II/base/task_result.h>

#include "../tests.h"



int
main()
{
  initlog();

  {
    // Create a background task that returns a value. Associate this
    // background task with a TaskResult object.
    Threads::TaskResult<int> t = Threads::new_task([]() {
      std::this_thread::sleep_for(std::chrono::seconds(1));
      return 42;
    });

    // Wait for that background task to finish, then get its
    // value:
    deallog << t.value() << std::endl;

    // Now set things to a different value, and output it:
    t.emplace_object(43);
    deallog << t.value() << std::endl;
  }

  {
    // Create an empty TaskResult object. Then set its value to
    // something and check what that is:
    Threads::TaskResult<int> t;
    t.emplace_object(44);
    deallog << t.value() << std::endl;
  }
}
