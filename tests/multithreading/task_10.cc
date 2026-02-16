// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2009 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// verify Threads::Task::joinable()

#include <deal.II/base/thread_management.h>

#include "../tests.h"


void
test()
{
  std::this_thread::sleep_for(std::chrono::seconds(1));
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
