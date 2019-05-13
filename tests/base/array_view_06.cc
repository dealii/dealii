// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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


// test that an ArrayView to an empty array view can actually be copied

#include <deal.II/base/array_view.h>

#include "../tests.h"


void
test()
{
  std::vector<int> v(10);

  ArrayView<int> a(&v[4], 0);
  ArrayView<int> b = a;

  ArrayView<int> c(nullptr, 0);
  ArrayView<int> d = a;

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
