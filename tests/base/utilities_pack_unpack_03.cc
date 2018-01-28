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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// test Utilities::pack/unpack on some types. this test checks that
// for trivially-copyable (small) types, packing is just a memcpy
// operation

#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/base/point.h>
#include <boost/serialization/utility.hpp>

#include <tuple>
#include <cstring>

struct X
{
  int i;
  int k;
  double d;

  bool operator == (const X &x) const
  {
    return i==x.i && k==x.k && d==x.d;
  }

  template <class Archive>
  void serialize(Archive &ar, const unsigned int)
  {
    ar   &i &k &d;
  }
};



template <typename T>
void check (const T &object)
{
  const std::vector<char> buffer = Utilities::pack (object);
  if (!(buffer.size() == sizeof(object))
      ||
      !(std::memcmp(buffer.data(), &object, buffer.size()) == 0)
      ||
      !(Utilities::unpack<T>(buffer) == object))
    deallog << "Fail!" << std::endl;
}


void test()
{
  check (std::make_pair(1, 3.14));
  check (X { 1, 2, 3.1415926 });
  check (std::tuple<int,double,char> {1,1,1});

  deallog << "OK!" << std::endl;
}

int main()
{
  initlog();

  test();
}
