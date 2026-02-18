// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test for class ArrayView

#include <deal.II/base/array_view.h>

#include "../tests.h"


void
test()
{
  {
    std::vector<int> arr = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    ArrayView<int>   view(arr); // writable view
    for (auto &el : view)
      ++el;

    int i = 0;
    for (auto &&it = arr.cbegin(); it != arr.cend(); ++it, ++i)
      AssertThrow(*it == i + 1, ExcInternalError());
  }

  {
    std::vector<int>     arr = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    const ArrayView<int> view(arr); // writable view
    for (auto &el : view)
      ++el;

    int i = 0;
    for (auto &&it = arr.cbegin(); it != arr.cend(); ++it, ++i)
      AssertThrow(*it == i + 1, ExcInternalError());
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
