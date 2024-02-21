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


// test the ArrayView constructor that converts from std::vector

#include <deal.II/base/array_view.h>

#include "../tests.h"


void
test()
{
  // converting a non-const vector to an ArrayView to const or
  // non-const data should work
  {
    std::vector<double>     v(10);
    ArrayView<double>       a1(v);
    ArrayView<const double> a2(v);
  }

  // converting a const vector to an ArrayView to const
  // data should work.
  //
  // converting to an ArrayView<double> will not work
  {
    const std::vector<double> v(10);
    ArrayView<const double>   a2(v);
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
