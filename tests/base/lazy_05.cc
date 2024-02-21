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
// Test move semantics of Lazy<T>.
//

#include <deal.II/base/lazy.h>

#include <iostream>

#include "../tests.h"


int
main()
{
  initlog();

  {
    Lazy<int> lazy_integer;
    lazy_integer.ensure_initialized([&]() {
      deallog << "... [ initialized ] ... " << std::endl;
      return 42;
    });
    deallog << "Object: " << lazy_integer.value() << std::endl;

    // Move constructor:
    Lazy<int> lazy_integer_2(std::move(lazy_integer));
    deallog << "Move result: " << lazy_integer_2.value() << std::endl;
    Assert(lazy_integer.has_value() == false, ExcInternalError());
  }

  {
    Lazy<int> lazy_integer;
    lazy_integer.ensure_initialized([&]() {
      deallog << "... [ initialized ] ... " << std::endl;
      return 42;
    });
    deallog << "Object: " << lazy_integer.value() << std::endl;

    // Move assignment
    Lazy<int> lazy_integer_2;
    lazy_integer_2 = std::move(lazy_integer);
    deallog << "Move result: " << lazy_integer_2.value() << std::endl;
    Assert(lazy_integer.has_value() == false, ExcInternalError());
  }


  deallog << "OK!" << std::endl;
  return 0;
}
