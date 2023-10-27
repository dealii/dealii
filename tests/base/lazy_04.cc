// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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

//
// Use Lazy<T> concurrently
//

#include <deal.II/base/lazy.h>

#include <iostream>

#include "../tests.h"

using namespace dealii;

int
main()
{
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
