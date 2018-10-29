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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Adoptation of base/utilities_pack_unpack_01.cc as a quick test.

// test Utilities::pack/unpack on some types.

#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <iostream>

using namespace dealii;

template <int dim>
void
test(const unsigned int &size)
{
  std::vector<Point<dim>> points(size);

  unsigned int i = 0;
  for (auto &p : points)
    for (unsigned int d = 0; d < dim; ++d)
      p[d] = 11 * d + size * 10 + 3 * i++;

  auto buffer = Utilities::pack(points);

  auto unpacked = Utilities::unpack<std::vector<Point<dim>>>(buffer);

  i       = 0;
  bool ok = true;
  for (const auto &p : points)
    if (p.distance(unpacked[i++]) > 1e-12)
      {
        std::cout << "NOT OK: " << p << " != " << unpacked[i - 1] << std::endl;
        ok = false;
      }

  AssertThrow(ok, ExcInternalError());
}

int
main()
{
  test<1>(10);
  test<2>(10);
  test<3>(10);
}
