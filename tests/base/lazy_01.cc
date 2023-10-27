// ---------------------------------------------------------------------
//
// Copyright (C) 2023 - 2023 by the deal.II authors
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
// Test the ensure_initialized() and value_or_initialize() interface of
// Lazy<T>.
//

#include <deal.II/base/lazy.h>

#include <iostream>

#include "../tests.h"

using namespace dealii;

int
main()
{
  initlog();

  {
    Lazy<int> lazy_integer;
    deallog << "lazy_integer.has_value() = " << lazy_integer.has_value()
            << std::endl;

    lazy_integer.ensure_initialized([&]() {
      deallog << "[initializing object]" << std::endl;
      return 42;
    });
    deallog << "lazy_integer.has_value() = " << lazy_integer.has_value()
            << std::endl;
    deallog << "lazy_integer.value() = " << lazy_integer.value() << std::endl;

    lazy_integer.ensure_initialized([&]() {
      deallog << "[initializing object] again... not good." << std::endl;
      return -1;
    });

    deallog << "lazy_integer.value() = " << lazy_integer.value() << std::endl;
  }

  {
    Lazy<int> lazy_integer;

    deallog << "lazy_integer.value() = "
            << lazy_integer.value_or_initialize([&]() {
                 deallog << "... [initializing object] ... ";
                 return 42;
               })
            << std::endl;
  }

  deallog << "OK!" << std::endl;
  return 0;
}
