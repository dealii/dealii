// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for AlignedVector<unsigned int> which tests the basic stuff in the
// aligned vector

#include <deal.II/base/aligned_vector.h>

#include <deque>

#include "../tests.h"

void
test()
{
  using VEC = AlignedVector<unsigned int>;
  VEC a(4);
  deallog << "Constructor: ";
  for (unsigned int i = 0; i < a.size(); ++i)
    deallog << a[i] << ' ';
  deallog << std::endl;

  a[2] = 1;
  a.push_back(5);
  a.push_back(42);

  VEC b(a);
  b.push_back(27);
  a.insert_back(b.begin(), b.end());

  // Check the range constructor with and without conversion.
  {
    VEC c(b.begin(), b.end());
    AssertThrow(c == b, ExcInternalError());

    std::vector<int> temp(b.begin(), b.end());
    VEC              d(temp.begin(), temp.end());
    AssertThrow(c == d, ExcInternalError());
  }

  // Check the range constructor with a random-access iterator which is not
  // contiguous
  {
    // Use a large enough deque to guarantee that we have multiple blocks
    std::deque<unsigned int> temp(8192);
    std::iota(temp.begin(), temp.end(), 0u);

    VEC e(temp.begin(), temp.end());
    AssertThrow(std::equal(e.begin(), e.end(), temp.begin()),
                ExcInternalError());
  }

  // Also check large insertions for equality and iterator position
  {
    std::deque<unsigned int> temp(8192);
    std::iota(temp.begin(), temp.end(), 0u);

    VEC        f(temp.begin(), temp.begin() + temp.size() / 4);
    const auto it0 =
      f.insert(f.end(), temp.begin() + temp.size() / 2, temp.end());
    AssertThrow(static_cast<std::size_t>(it0 - f.begin()) == temp.size() / 4,
                ExcInternalError());
    AssertThrow(*it0 == temp[temp.size() / 2], ExcInternalError());
    AssertThrow(*(it0 - 1) == temp[temp.size() / 4 - 1], ExcInternalError());
    AssertThrow(f.back() == temp.back(), ExcInternalError());

    const auto it1 = f.insert(f.begin() + temp.size() / 4,
                              temp.begin() + temp.size() / 4,
                              temp.begin() + temp.size() / 2);
    AssertThrow(static_cast<std::size_t>(it1 - f.begin()) == temp.size() / 4,
                ExcInternalError());
    AssertThrow(*it1 == temp[temp.size() / 4], ExcInternalError());
    AssertThrow(std::equal(f.begin(), f.end(), temp.begin()),
                ExcInternalError());
  }

  deallog << "back Insertion: ";
  for (unsigned int i = 0; i < a.size(); ++i)
    deallog << a[i] << ' ';
  deallog << std::endl;

  {
    AlignedVector<unsigned short> temp(4);
    std::fill(temp.begin(), temp.end(), 42u);
    a.insert(a.begin() + 4, temp.begin(), temp.end());
  }
  deallog << "Insertion at position 4: ";

  for (unsigned int i = 0; i < a.size(); ++i)
    deallog << a[i] << ' ';
  deallog << std::endl;

  deallog << "Memory Shrinking: ";
  deallog << a.memory_consumption() << " to ";
  a.resize(4);
  a.shrink_to_fit();
  deallog << a.memory_consumption() << std::endl;
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
  a.resize(100000, 1);
  deallog << "Check large initialization: ";
  for (unsigned int i = 0; i < 100000; ++i)
    AssertDimension(a[i], 1);
  deallog << "OK" << std::endl;

  // check resize for large vectors
  deallog << "Check large resize: ";
  a.resize(200000, 2);
  a.resize(400000);
  for (unsigned int i = 0; i < 100000; ++i)
    AssertDimension(a[i], 1);
  for (unsigned int i = 100000; i < 200000; ++i)
    AssertDimension(a[i], 2);
  for (unsigned int i = 200000; i < 400000; ++i)
    AssertDimension(a[i], 0);
  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
