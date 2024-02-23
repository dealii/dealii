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


// testing serialize function for class BoundingBox

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include "../tests.h"

template <int spacedim>
void
test(const unsigned int &size)
{
  std::vector<BoundingBox<spacedim>> b_boxes(size);

  for (auto &b : b_boxes)
    {
      auto p1 = random_point<spacedim>();
      auto p2 = p1;
      for (unsigned int i = 0; i < spacedim; ++i)
        p2[i] += i;
      b = BoundingBox<spacedim>(std::make_pair(p1, p2));
    }

  auto buffer = Utilities::pack(b_boxes);

  auto unpacked = Utilities::unpack<std::vector<BoundingBox<spacedim>>>(buffer);

  unsigned int i  = 0;
  bool         ok = true;
  for (auto &b : b_boxes)
    {
      const auto &b_points = b.get_boundary_points();
      const auto &u_points = unpacked[i++].get_boundary_points();
      if (b_points.first.distance(u_points.first) > 1e-12)
        {
          deallog << "NOT OK: " << b_points.first << " != " << u_points.first
                  << std::endl;
          ok = false;
        }
      else if (b_points.second.distance(u_points.second) > 1e-12)
        {
          deallog << "NOT OK: " << b_points.second << " != " << u_points.second
                  << std::endl;
          ok = false;
        }
    }

  if (ok)
    deallog << "OK!" << std::endl;
}

int
main()
{
  initlog();

  test<1>(10);
  test<2>(20);
  test<3>(30);
}
