// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
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
  std::vector<int> v(10);

  ArrayView<int> a(&v[4], 3); // writable view
  a[2] = 42;

  Assert(a[2] == 42, ExcInternalError());
  Assert(v[6] == 42, ExcInternalError());

  ArrayView<const int> a2(&v[4], 3); // readable view
  Assert(a2[2] == 42, ExcInternalError());

  ArrayView<const int> a3(a); // readable view, converted from 'a'
  Assert(a3[2] == 42, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
