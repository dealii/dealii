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
// Use Lazy<T> in const context:
//

#include <deal.II/base/lazy.h>

#include <iostream>

#include "../tests.h"


class Container
{
public:
  int
  get() const
  {
    lazy_integer.ensure_initialized([]() { return 42; });

    if (lazy_integer.has_value())
      return lazy_integer.value();

    // unreachable
    return lazy_integer.value_or_initialize([]() { return 43; });
  }

private:
  Lazy<int> lazy_integer;
};

int
main()
{
  initlog();

  Container foo;

  foo.get();

  deallog << "OK!" << std::endl;
  return 0;
}
