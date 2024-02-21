// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
