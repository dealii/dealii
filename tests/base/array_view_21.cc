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


// test for class ArrayView initialized from a std::initializer_list.

#include <deal.II/base/array_view.h>

#include "../tests.h"


void
test(const ArrayView<const int> &a)
{
  deallog << "Size=" << a.size() << std::endl;
  for (const auto &i : a)
    deallog << "  " << i << std::endl;
}


int
main()
{
  initlog();

  test({});
  test({1, 2, 3});
}
