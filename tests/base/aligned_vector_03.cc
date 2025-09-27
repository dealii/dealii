// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for AlignedVector<bool>

#include <deal.II/base/aligned_vector.h>

#include "../tests.h"


void
test()
{
  using VEC = AlignedVector<bool>;
  VEC a(4);
  deallog << "Constructor: ";
  for (unsigned int i = 0; i < a.size(); ++i)
    deallog << a[i] << ' ';
  deallog << std::endl;

  a[2] = true;
  a.push_back(true);
  a.push_back(false);

  VEC b(a);
  b.push_back(true);
  a.insert_back(b.begin(), b.end());

  deallog << "Insertion: ";
  for (unsigned int i = 0; i < a.size(); ++i)
    deallog << a[i] << ' ';
  deallog << std::endl;

  a.resize(4);
  deallog << "Shrinking: ";
  for (unsigned int i = 0; i < a.size(); ++i)
    deallog << a[i] << ' ';
  deallog << std::endl;

  a.reserve(100);
  deallog << "Reserve: ";
  for (unsigned int i = 0; i < a.size(); ++i)
    deallog << a[i] << ' ';
  deallog << std::endl;

  a = b;
  deallog << "Assignment: ";
  for (unsigned int i = 0; i < a.size(); ++i)
    deallog << a[i] << ' ';
  deallog << std::endl;

  // check setting elements for large vectors
  a.resize(0);
  a.resize(100000, true);
  deallog << "Check large initialization: ";
  for (unsigned int i = 0; i < 100000; ++i)
    AssertDimension(a[i], true);
  deallog << "OK" << std::endl;

  // check resize for large vectors
  deallog << "Check large resize: ";
  a.resize(200000, false);
  a.resize(400000, true);
  for (unsigned int i = 0; i < 100000; ++i)
    AssertDimension(a[i], true);
  for (unsigned int i = 100000; i < 200000; ++i)
    AssertDimension(a[i], false);
  for (unsigned int i = 200000; i < 400000; ++i)
    AssertDimension(a[i], true);
  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
