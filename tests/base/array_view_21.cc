// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
