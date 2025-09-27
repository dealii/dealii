// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

//
// Use Lazy<T> concurrently
//

#include <deal.II/base/lazy.h>

#include <iostream>

#include "../tests.h"


int
main()
{
  return 1; // test is disabled for now as it hangs randomly
  initlog();

  Lazy<int> lazy_integer;

  const auto foo = [&]() {
    lazy_integer.ensure_initialized([&]() {
      deallog << "... [ initialized ] ... " << std::endl;
      return 42;
    });
  };

  Threads::Task<> t1 = Threads::new_task(foo);
  Threads::Task<> t2 = Threads::new_task(foo);
  Threads::Task<> t3 = Threads::new_task(foo);
  Threads::Task<> t4 = Threads::new_task(foo);
  t4.join();
  t3.join();
  t2.join();
  t1.join();

  deallog << "Object: " << lazy_integer.value() << std::endl;

  deallog << "OK!" << std::endl;
  return 0;
}
