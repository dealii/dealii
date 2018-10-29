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


// test Utilities::pack/unpack on some types.

#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include "../tests.h"

template <int dim>
void
test(const unsigned int &size)
{
  std::vector<Point<dim>> points(size);

  for (auto &p : points)
    p = random_point<dim>();

  auto buffer = Utilities::pack(points);

  auto unpacked = Utilities::unpack<std::vector<Point<dim>>>(buffer);

  unsigned int i  = 0;
  bool         ok = true;
  for (const auto &p : points)
    if (p.distance(unpacked[i++]) > 1e-12)
      {
        deallog << "NOT OK: " << p << " != " << unpacked[i - 1] << std::endl;
        ok = false;
      }

  if (ok)
    deallog << "OK!" << std::endl;
}

int
main()
{
  initlog();

  test<1>(10);
  test<2>(10);
  test<3>(10);
}
