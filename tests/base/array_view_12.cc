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
