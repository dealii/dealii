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


// verify that TaskGroup::size() does what we want

#include <deal.II/base/thread_management.h>

#include "../tests.h"


void
test(int i)
{
  std::this_thread::sleep_for(std::chrono::seconds(1));
}



int
main()
{
  initlog();

  Threads::TaskGroup<> tg;
  tg += Threads::new_task(test, 1);
  deallog << tg.size() << std::endl;

  tg += Threads::new_task(test, 2);
  deallog << tg.size() << std::endl;

  tg.join_all();
  deallog << tg.size() << std::endl;
}
