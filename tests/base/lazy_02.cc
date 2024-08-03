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
// Test copy, move, and reset semantics of Lazy<T>.
//

#include <deal.II/base/lazy.h>

#include <iostream>

#include "../tests.h"


int
main()
{
  initlog();

  Lazy<int> lazy_integer;
  lazy_integer.ensure_initialized([&]() {
    deallog << "... [ initialized ] ... " << std::endl;
    return 42;
  });
  deallog << "Object: " << lazy_integer.value() << std::endl;

  Lazy<int> lazy_integer_2(lazy_integer);
  deallog << "Copied: " << lazy_integer_2.value() << std::endl;

  Lazy<int> lazy_integer_3;
  lazy_integer_3 = lazy_integer_2;
  deallog << "Assigned: " << lazy_integer_3.value() << std::endl;

  Lazy<int> lazy_integer_4 = [&]() {
    Lazy<int> temp = lazy_integer;
    return temp;
  }();
  deallog << "Moved: " << lazy_integer_4.value() << std::endl;

  lazy_integer_4 = [&]() {
    Lazy<int> temp = lazy_integer;
    return temp;
  }();
  deallog << "Moved assigned: " << lazy_integer_4.value() << std::endl;

  deallog << "... [ reset ] ... " << std::endl;
  lazy_integer_4.reset();
  lazy_integer_4.ensure_initialized([&]() {
    deallog << "... [ initialized ] ... " << std::endl;
    return 43;
  });

  deallog << "... [ reset ] ... " << std::endl;
  lazy_integer_4.reset();
  lazy_integer_4.ensure_initialized([&]() {
    deallog << "... [ initialized ] ... " << std::endl;
    return 43;
  });
  deallog << "Reset: " << lazy_integer_4.value() << std::endl;

  deallog << "... [ reset ] ... " << std::endl;
  lazy_integer_4.reset();
  deallog << "value_or_initialize: " +
               std::to_string(lazy_integer_4.value_or_initialize([&]() {
                 deallog << "... [ initialized ] ... " << std::endl;
                 return 45;
               }))
          << std::endl;

  deallog << "OK!" << std::endl;
  return 0;
}
