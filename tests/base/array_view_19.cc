// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for class ArrayView with AlignedVector, using make_array_view

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/array_view.h>

#include "../tests.h"


void
test()
{
  AlignedVector<int> v(10, 1.0);

  {
    ArrayView<int> a(&v[4], 3); // writable view
    a[2] = 42;

    Assert(a[2] == 42, ExcInternalError());
    Assert(v[6] == 42, ExcInternalError());
  }

  {
    ArrayView<int> a = make_array_view(v); // writable view
    a[1]             = 43;

    Assert(v[1] == 43, ExcInternalError());
    Assert(a[6] == 42, ExcInternalError());
  }

  {
    ArrayView<const int> a2(&v[4], 3); // readable view
    Assert(a2[2] == 42, ExcInternalError());
  }

  {
    ArrayView<const int> a2 = make_array_view(v); // readable view
    Assert(a2[1] == 43, ExcInternalError());
    Assert(a2[6] == 42, ExcInternalError());
  }

  {
    // writable view on subrange
    ArrayView<int> a2 = make_array_view(v, 1, 7);
    a2[0]             = 44;
    Assert(a2[0] == 44, ExcInternalError());
    Assert(a2[5] == 42, ExcInternalError());
  }

  {
    // readable view on subrange
    ArrayView<const int> a2 = make_array_view(v, 1, 7);
    Assert(a2[0] == 44, ExcInternalError());
    Assert(a2[5] == 42, ExcInternalError());
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
