// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Extract a representative vector of BoundingBox objects from an RTree

#include <deal.II/base/patterns.h>

#include <deal.II/boost_adaptors/bounding_box.h>
#include <deal.II/boost_adaptors/point.h>

#include <deal.II/numerics/rtree.h>

#include <algorithm>

#include "../tests.h"

int
main()
{
  initlog(true);

  const unsigned int    N = 20;
  std::vector<Point<2>> points(N);

  for (auto &p : points)
    p = random_point<2>();

  // Make sure we have at most two points in each box
  const auto tree = pack_rtree<boost::geometry::index::linear<2>>(points);

  const auto boxes = extract_rtree_level(tree, 1);
  deallog << "LEVEL 1:  N boxes: " << boxes.size() << std::endl;

  for (const auto &b : boxes)
    deallog << "Box: " << Patterns::Tools::to_string(b.get_boundary_points())
            << std::endl;
}
