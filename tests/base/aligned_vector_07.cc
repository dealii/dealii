// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2021 by the deal.II authors
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
