// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test Utilities::pack/unpack on some types. this test checks that
// for trivially-copyable (small) types, packing is just a memcpy
// operation

#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <boost/serialization/utility.hpp>

#include <array>

#include "../tests.h"

struct X
{
  int    i;
  int    k;
  double d;

  bool
  operator!=(const X &x) const
  {
    return i != x.i || k != x.k || d != x.d;
  }

  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int)
  {
    ar &i &k &d;
  }
};



template <typename T>
void
check(const T &object)
{
  const std::vector<char> buffer = Utilities::pack(object);
  if (buffer.size() != sizeof(object))
    deallog << buffer.size() << " should be " << sizeof(object) << '!'
            << std::endl;
  else
    deallog << "same size!" << std::endl;

  if (std::memcmp(buffer.data(), &object, buffer.size()) != 0)
    deallog << "std::memcmp failed!" << std::endl;
  else
    deallog << "std::memcmp passed!" << std::endl;

  if (Utilities::unpack<T>(buffer) != object)
    deallog << "Comparing the objects failed!" << std::endl;
  else
    deallog << "Comparing the objects passed!" << std::endl;

  deallog << std::endl;
}


void
test()
{
  deallog << "std::array:" << std::endl;
  check(std::array<int, 3>{{1, 2, 3}});
  deallog << "struct X:" << std::endl;
  check(X{1, 2, numbers::PI});
  deallog << "double:" << std::endl;
  check(1.);
}

int
main()
{
  initlog();

  test();
}
