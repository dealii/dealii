// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Adaptation of base/utilities_pack_unpack_01.cc as a quick test.

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
