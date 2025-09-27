// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for class ArrayView. check make_array_view for Table arguments

#include <deal.II/base/array_view.h>
#include <deal.II/base/table.h>

#include "../tests.h"


void
test()
{
  Table<2, int> v(8, 10);

  ArrayView<int> a =
    make_array_view(v, 4, 7, 3); // writable view to part of a row
  a[2] = 42;

  Assert(a[2] == 42, ExcInternalError());
  Assert(v[4][9] == 42, ExcInternalError());

  ArrayView<const int> a2 = make_array_view(v, 4, 7, 3); // readable view
  Assert(a2[2] == 42, ExcInternalError());


  // also check a different way of creating a readable view
  ArrayView<const int> a3 =
    make_array_view(const_cast<const Table<2, int> &>(v), 4, 7, 3);
  Assert(a3[2] == 42, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
