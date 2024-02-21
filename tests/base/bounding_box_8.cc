// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for BoundingBox<int dim> which tests the function
// has_overlap_with()

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/point.h>

#include "../tests.h"

BoundingBox<1>
generate_bbox(const double left, const double right)
{
  Point<1> p1;
  Point<1> p2;

  p1[0] = left;
  p2[0] = right;

  BoundingBox<1> bbox(std::make_pair(p1, p2));

  return bbox;
}

void
test_bounding_box(const double left_a,
                  const double right_a,
                  const double left_b,
                  const double right_b,
                  const double tolerance)
{
  const auto bbox_a = generate_bbox(left_a, right_a);
  const auto bbox_b = generate_bbox(left_b, right_b);

  deallog << "Bounding box A: ";
  deallog << "[" << bbox_a.get_boundary_points().first << ", "
          << bbox_a.get_boundary_points().second << "]" << std::endl;
  deallog << "Bounding box B: ";
  deallog << "[" << bbox_b.get_boundary_points().first << ", "
          << bbox_b.get_boundary_points().second << "]" << std::endl;
  deallog << "Has overlap with: " << bbox_a.has_overlap_with(bbox_b, tolerance)
          << std::endl;
}

int
main()
{
  initlog();

  test_bounding_box(1.0, 2.0, 2.0, 3.0, 1e-12);
  test_bounding_box(1.0, 2.0, 2.0 + 1e-11, 3.0, 1e-12);
  test_bounding_box(1.0, 2.0, 2.0 + 1e-11, 3.0, 1e-10);

  test_bounding_box(-1.0, 0.0, 0.0, 1.0, 1e-12);
  test_bounding_box(-1.0, 0.0, 0.0 + 1e-11, 1.0, 1e-12);
  test_bounding_box(-1.0, 0.0, 0.0 + 1e-11, 1.0, 1e-10);

  test_bounding_box(-2.0, -1.0, -1.0, 0.0, 1e-12);
  test_bounding_box(-2.0, -1.0, -1.0 + 1e-11, 0.0, 1e-12);
  test_bounding_box(-2.0, -1.0, -1.0 + 1e-11, 0.0, 1e-10);

  return 0;
}
