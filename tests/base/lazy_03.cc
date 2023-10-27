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
// Use Lazy<T> in const context:
//

#include <deal.II/base/lazy.h>

#include <iostream>

#include "../tests.h"

using namespace dealii;

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
