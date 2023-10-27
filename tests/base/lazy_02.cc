// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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
// Test copy, move, and reset semantics of Lazy<T>.
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
  deallog << "value_or_initialize: "
          << lazy_integer_4.value_or_initialize([&]() {
               deallog << "... [ initialized ] ... ";
               return 45;
             })
          << std::endl;

  deallog << "OK!" << std::endl;
  return 0;
}
