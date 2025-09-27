// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test for AlignedVector<AlignedVector<unsigned int> >, which
// uncovered a difficult-to-track bug in the move assignment operator
// related to the destruction of element objects.

#include <deal.II/base/aligned_vector.h>

#include "../tests.h"



void
test()
{
  AlignedVector<unsigned int> a(1);
  a[0] = 42;

  AlignedVector<AlignedVector<unsigned int>> v(1);
  v[0] = a;

  deallog << "Resizing..." << std::endl;
  v.resize(2);
  a[0] = 932;
  v[1] = a;

  deallog << "Shrinking..." << std::endl;
  v.resize(1);

  deallog << "Checking: size=" << v[0].size() << ", element=" << v[0][0]
          << std::endl;

  deallog << "Clearing..." << std::endl;
  v.clear();

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
