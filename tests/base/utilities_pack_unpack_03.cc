// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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


// test Utilities::pack/unpack on some types. this test checks that
// for trivially-copyable (small) types, packing is just a memcpy
// operation

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
    deallog << buffer.size() << " should be " << sizeof(object) << "!"
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
  check(X{1, 2, 3.1415926});
  deallog << "double:" << std::endl;
  check(1.);
}

int
main()
{
  initlog();

  test();
}
