// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// test for class ArrayView. check make_array_view for whole vectors

#include <deal.II/base/array_view.h>

#include "../tests.h"


void
test()
{
  std::vector<int> v(10);

  ArrayView<int> a = make_array_view(v); // writable view
  a[2]             = 42;

  Assert(a[2] == 42, ExcInternalError());
  Assert(v[2] == 42, ExcInternalError());

  ArrayView<const int> a2 = make_array_view(v); // readable view
  Assert(a2[2] == 42, ExcInternalError());

  ArrayView<const int> a3 =
    make_array_view(const_cast<const std::vector<int> &>(v));
  Assert(a3[2] == 42, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
