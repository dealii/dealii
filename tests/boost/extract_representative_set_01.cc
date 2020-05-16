// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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
  initlog(1);

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