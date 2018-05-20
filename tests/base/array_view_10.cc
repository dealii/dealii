// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

// test the ArrayView constructor that converts from std::vector

#include "../tests.h"

#include <deal.II/base/array_view.h>

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
