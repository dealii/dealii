// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for class ArrayView

#include <deal.II/base/array_view.h>

#include "../tests.h"


void
test()
{
  std::vector<int>     arr1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  ArrayView<int>       view1(arr1);
  ArrayView<const int> view1_const(arr1);

  std::vector<int>     arr2 = {0, 1, 2, 3, 4};
  ArrayView<int>       view2(arr2);
  ArrayView<const int> view2_const(arr2);

  deallog << std::boolalpha;

  deallog << "view1==view1: " << (view1 == view1) << '\n';
  deallog << "view1!=view1: " << (view1 != view1) << '\n';
  deallog << "view1==view1_const: " << (view1 == view1_const) << '\n';
  deallog << "view1!=view1_const: " << (view1 != view1_const) << '\n';
  deallog << "view1==view2: " << (view1 == view2) << '\n';
  deallog << "view1!=view2: " << (view1 != view2) << '\n';
  deallog << "view1==view2_const: " << (view1 == view2_const) << '\n';
  deallog << "view1!=view2_const: " << (view1 != view2_const) << std::endl;

  deallog << "view1_const==view1: " << (view1_const == view1) << '\n';
  deallog << "view1_const!=view1: " << (view1_const != view1) << '\n';
  deallog << "view1_const==view1_const: " << (view1_const == view1_const)
          << '\n';
  deallog << "view1_const!=view1_const: " << (view1_const != view1_const)
          << '\n';
  deallog << "view1_const==view2: " << (view1_const == view2) << '\n';
  deallog << "view1_const!=view2: " << (view1_const != view2) << '\n';
  deallog << "view1_const==view2_const: " << (view1_const == view2_const)
          << '\n';
  deallog << "view1_const!=view2_const: " << (view1_const != view2_const)
          << std::endl;

  deallog << "view2==view1: " << (view2 == view1) << '\n';
  deallog << "view2!=view1: " << (view2 != view1) << '\n';
  deallog << "view2==view1_const: " << (view2 == view1_const) << '\n';
  deallog << "view2!=view1_const: " << (view2 != view1_const) << '\n';
  deallog << "view2==view2: " << (view2 == view2) << '\n';
  deallog << "view2!=view2: " << (view2 != view2) << '\n';
  deallog << "view2==view2_const: " << (view2 == view2_const) << '\n';
  deallog << "view2!=view2_const: " << (view2 != view2_const) << std::endl;

  deallog << "view2_const==view1: " << (view2_const == view1) << '\n';
  deallog << "view2_const!=view1: " << (view2_const != view1) << '\n';
  deallog << "view2_const==view1_const: " << (view2_const == view1_const)
          << '\n';
  deallog << "view2_const!=view1_const: " << (view2_const != view1_const)
          << '\n';
  deallog << "view2_const==view2: " << (view2_const == view2) << '\n';
  deallog << "view2_const!=view2: " << (view2_const != view2) << '\n';
  deallog << "view2_const==view2_const: " << (view2_const == view2_const)
          << '\n';
  deallog << "view2_const!=view2_const: " << (view2_const != view2_const)
          << std::endl;
}



int
main()
{
  initlog();

  test();
}
