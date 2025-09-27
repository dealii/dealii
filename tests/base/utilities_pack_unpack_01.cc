// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
