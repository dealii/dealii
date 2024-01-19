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
