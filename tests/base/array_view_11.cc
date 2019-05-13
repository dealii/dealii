// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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
