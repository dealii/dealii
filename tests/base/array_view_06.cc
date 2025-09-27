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
