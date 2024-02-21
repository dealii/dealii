// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for BoundingBox<unsigned int spacedim> which test basic stuff:
// creation in various dimensions and checking for points inside

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/point.h>

#include "../tests.h"

template <int spacedim>
void
test_bounding_box()
{
  BoundingBox<spacedim> a;
  deallog << "Empty constructor: " << std::endl;
  deallog << a.get_boundary_points().first << std::endl;
  deallog << a.get_boundary_points().second << std::endl;

  std::pair<Point<spacedim>, Point<spacedim>> boundaries;

  for (int i = 0; i < spacedim; ++i)
    {
      boundaries.first[i]  = 0.2 - i * 0.2;
      boundaries.second[i] = 0.8 + i * 0.8;
    }

  BoundingBox<spacedim> b(boundaries);
  deallog << "Boundary points: " << std::endl;
  deallog << b.get_boundary_points().first << std::endl;
  deallog << b.get_boundary_points().second << std::endl;

  deallog << "Boundary points are inside: " << b.point_inside(boundaries.first)
          << ' ' << b.point_inside(boundaries.second) << std::endl;

  std::vector<Point<spacedim>> test_points;

  // To guarantee points are inside we take a convex combination
  double c;
  for (int t = 0; t < 5; ++t)
    {
      Point<spacedim> test_pt;
      c = (1 + t) / 5.0;
      for (unsigned int i = 0; i < spacedim; ++i)
        test_pt[i] = boundaries.first[i] * c + boundaries.second[i] * (1 - c);

      test_points.push_back(test_pt);
    }

  deallog << "Points inside: " << std::endl;
  for (unsigned int i = 0; i < test_points.size(); ++i)
    deallog << test_points[i]
            << " is inside: " << b.point_inside(test_points[i]) << std::endl;

  deallog << std::endl;

  // Verify that we get the same box when we use all those points together
  {
    test_points.push_back(boundaries.first);
    test_points.push_back(boundaries.second);

    BoundingBox<spacedim> b2(test_points);
    deallog << "Boxes should be equal : " << (b2 == b) << std::endl;
  }

  test_points.clear();

  // To create outside points we take a non-convex combination
  for (int t = 0; t < 5; ++t)
    {
      Point<spacedim> test_pt;
      c = (1 + t) * 2.5;
      if (t % 2 == 0)
        c = -c; // Changing the sign sometimes..
      for (unsigned int i = 0; i < spacedim; ++i)
        test_pt[i] = boundaries.first[i] * c + boundaries.second[i] * c;

      test_points.push_back(test_pt);
    }

  deallog << "Points outside:" << std::endl;
  for (unsigned int i = 0; i < test_points.size(); ++i)
    deallog << test_points[i]
            << " is inside: " << b.point_inside(test_points[i]) << std::endl;

  // Similarly, verify that we get different boxes since some points are
  // outside:
  {
    BoundingBox<spacedim> b2(test_points);
    deallog << "Boxes should not be equal : " << (b2 != b) << std::endl;
  }
  deallog << std::endl;

  // Initialize box with point
  {
    Point<spacedim> p;
    for (unsigned int i = 0; i < spacedim; ++i)
      p[i] = i + 1;

    BoundingBox<spacedim> b(p);
    deallog << "Boundary points: " << std::endl;
    deallog << b.get_boundary_points().first << std::endl;
    deallog << b.get_boundary_points().second << std::endl;
  }
  deallog << std::endl;

  // Initialize box with box
  {
    BoundingBox<spacedim> bb(b);
    deallog << "Boundary points: " << std::endl;
    deallog << bb.get_boundary_points().first << std::endl;
    deallog << bb.get_boundary_points().second << std::endl;
  }
  deallog << std::endl;

  // Initialize box with box
  {
    BoundingBox<spacedim> bb;
    bb = b;
    deallog << "Boundary points: " << std::endl;
    deallog << bb.get_boundary_points().first << std::endl;
    deallog << bb.get_boundary_points().second << std::endl;
  }
  deallog << std::endl;
}

void
test_unitary()
{
  std::pair<Point<3>, Point<3>> boundaries;

  for (int i = 0; i < 3; ++i)
    {
      boundaries.second[i] = 1.0;
    }

  BoundingBox<3> b(boundaries);
  deallog << "Boundary points:" << std::endl;
  deallog << b.get_boundary_points().first << std::endl;
  deallog << b.get_boundary_points().second << std::endl;

  Point<3> p1(1.0, 0, 0);
  Point<3> p2(0, 1.0, 0);
  Point<3> p3(0, 0, 1.0);

  deallog << "Checking if all vertices are inside: "
          << b.point_inside(boundaries.first) << ' '
          << b.point_inside(boundaries.second) << std::endl;

  deallog << b.point_inside(p1) << ' ' << b.point_inside(p2) << ' '
          << b.point_inside(p3) << ' ' << b.point_inside(p1 + p2) << ' '
          << b.point_inside(p2 + p3) << ' ' << b.point_inside(p1 + p3) << ' '
          << std::endl;

  double eps = std::numeric_limits<double>::epsilon();
  AssertThrow(b.point_inside(Point<3>(0.0, 0.0, 1.0 + 1.0 * eps), 10. * eps) ==
                true,
              ExcMessage("failed."));
  AssertThrow(b.point_inside(Point<3>(0.0, 0.0, 1.0 + 10. * eps), 1.0 * eps) ==
                false,
              ExcMessage("failed."));
  AssertThrow(b.point_inside(Point<3>(0.0 - 1.0 * eps, 0.0, 0.0), 10. * eps) ==
                true,
              ExcMessage("failed."));
  AssertThrow(b.point_inside(Point<3>(0.0 - 10. * eps, 0.0, 0.0), 1.0 * eps) ==
                false,
              ExcMessage("failed."));
}

int
main()
{
  initlog();

  deallog << "Test: Bounding Box class " << std::endl;
  deallog << "Test for dimension 1" << std::endl;
  test_bounding_box<1>();

  deallog << std::endl << "Test for dimension 2" << std::endl;
  test_bounding_box<2>();

  deallog << std::endl << "Test for dimension 3" << std::endl;
  test_bounding_box<3>();

  deallog << std::endl << "Test for dimension 3, unitary box" << std::endl;
  test_unitary();
}
