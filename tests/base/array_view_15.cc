// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// test for class ArrayView for scalar type

#include <deal.II/base/array_view.h>

#include "../tests.h"


void
test()
{
  int v = 9;

  ArrayView<int> a1(v); // writable view

  Assert(a1.size() == 1, ExcInternalError());
  Assert(a1[0] == 9, ExcInternalError());
  v = 10;
  Assert(a1[0] == 10, ExcInternalError());
  a1[0] = 11;
  Assert(a1[0] == 11, ExcInternalError());

  ArrayView<int> a2(v); // writable view

  Assert(a2.size() == 1, ExcInternalError());
  Assert(a2[0] == 11, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
