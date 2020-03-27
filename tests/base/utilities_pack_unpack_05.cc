// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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


// test Utilities::pack/unpack on a consecutively built buffer of
// different types, here an array of doubles and a dealii::Point object.
// (based upon "utilities_pack_unpack_04")


#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include "../tests.h"



template <int N, int dim>
void
check(const double (&array)[N], const Point<dim>(&point))
{
  std::vector<char> buffer;

  // PACK BUFFER
  // add first object to buffer and store buffer size for later separation
  const std::size_t buffer_separator = Utilities::pack(array, buffer);
  // add second object to buffer
  Utilities::pack(point, buffer);

  // UNPACK BUFFER
  double unpacked_array[N];
  Utilities::unpack(buffer.cbegin(),
                    buffer.cbegin() + buffer_separator,
                    unpacked_array);
  Point<dim> unpacked_point =
    Utilities::unpack<Point<dim>>(buffer.cbegin() + buffer_separator,
                                  buffer.cend());

  // TEST RESULTS
  bool equal_array = true;
  for (unsigned int i = 0; i < N; ++i)
    if (array[i] != unpacked_array[i])
      {
        equal_array = false;
        break;
      }
  deallog << "compare array: " << (equal_array ? "OK" : "Failed") << std::endl;

  deallog << "compare point: "
          << (point.distance(unpacked_point) < 1e-12 ? "OK" : "Failed")
          << std::endl;
}


void
test()
{
  // try small arrays that are packed by just using memcpy
  Point<3> p1    = random_point<3>();
  double   x1[3] = {1, 2, 3};
  check(x1, p1);

  // now try much larger arrays that will actually be serialized
  // using BOOST
  const unsigned int N  = 10000;
  Point<N>           p2 = random_point<N>();
  double             x2[N];
  for (unsigned int i = 0; i < N; ++i)
    x2[i] = i;
  check(x2, p2);

  deallog << "OK!" << std::endl;
}

int
main()
{
  initlog();

  test();
}
