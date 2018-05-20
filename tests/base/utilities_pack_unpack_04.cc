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

// test Utilities::pack/unpack on some types. specifically, make sure
// that it works on arrays of primitive types

#include "../tests.h"

#include <boost/serialization/utility.hpp>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <cstring>
#include <tuple>

template <typename T, int N>
void
check(const T (&object)[N])
{
  const std::vector<char> buffer = Utilities::pack(object);
  T                       unpacked[N];
  Utilities::unpack(buffer, unpacked);

  const bool equal_sizes = (buffer.size() == sizeof(T) * N);
  deallog << "Buffer size check: " << (equal_sizes ? "OK" : "Failed")
          << std::endl;
  if(equal_sizes)
    deallog << "memcmp check: "
            << (std::memcmp(buffer.data(), &object, buffer.size()) == 0 ?
                  "OK" :
                  "Failed")
            << std::endl;

  bool equal = true;
  for(unsigned int i = 0; i < N; ++i)
    if(unpacked[i] != object[i])
      {
        equal = false;
        break;
      }
  deallog << "direct cmp: " << (equal ? "OK" : "Failed") << std::endl;
}

void
test()
{
  // try a small array that is packed by just using memcpy
  double x[3] = {1, 2, 3};
  check(x);

  // now try a much larger array that will actually be serialized
  // using BOOST
  const unsigned int N = 10000;
  double             y[N];
  for(unsigned int i = 0; i < N; ++i)
    y[i] = i;
  check(y);

  deallog << "OK!" << std::endl;
}

int
main()
{
  initlog();

  test();
}
