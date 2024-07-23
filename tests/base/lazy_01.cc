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
// Test the ensure_initialized() and value_or_initialize() interface of
// Lazy<T>.
//

#include <deal.II/base/lazy.h>

#include <iostream>

#include "../tests.h"


int
main()
{
  initlog();

  {
    deallog << "Part 1:" << std::endl;

    Lazy<int> lazy_integer;
    deallog << "lazy_integer.has_value() = " << lazy_integer.has_value()
            << std::endl;

    lazy_integer.ensure_initialized([&]() {
      deallog << "[initializing object]" << std::endl;
      return 42;
    });

    // Lazy computes the result on a different task, which may be on a
    // different thread. To avoid producing unreliable output, make
    // sure we don't say
    //   deallog << "lazy_integer.has_value() = " << lazy_integer.has_value()
    // but the following, where the order of computation and output are
    // clear and well defined:
    deallog << "lazy_integer.has_value() = " +
                 std::to_string(lazy_integer.has_value())
            << std::endl;
    deallog << "lazy_integer.value() = " << lazy_integer.value() << std::endl;

    lazy_integer.ensure_initialized([&]() {
      deallog << "[initializing object] again... not good." << std::endl;
      return -1;
    });

    deallog << "lazy_integer.value() = " + std::to_string(lazy_integer.value())
            << std::endl;
  }

  {
    deallog << "Part 2:" << std::endl;

    Lazy<int> lazy_integer;

    deallog << "lazy_integer.value() = " +
                 std::to_string(lazy_integer.value_or_initialize([&]() {
                   deallog << "... [initializing object] ... " << std::endl;
                   return 42;
                 }))
            << std::endl;
  }

  deallog << "OK!" << std::endl;
  return 0;
}
