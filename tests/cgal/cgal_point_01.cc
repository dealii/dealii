// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test conversion from CGAL to deal.II Point and vice-versa

#include <deal.II/base/point.h>

#include <deal.II/cgal/utilities.h>

#include <CGAL/IO/io.h>
#include <CGAL/Simple_cartesian.h>

#include "../tests.h"

using namespace CGALWrappers;

int
main()
{
  initlog();
  using CGALPoint = CGAL::Point_3<CGAL::Simple_cartesian<double>>;
  // Test conversion from deal.II Point to CGAL Point
  {
    const Point<3> dealii_point(1.0, 2.0, 3.0);
    const auto cgal_point = dealii_point_to_cgal_point<CGALPoint>(dealii_point);
    deallog << "CGAL Point: " << cgal_point << std::endl;
  }

  // Test conversion from CGAL Point to deal.II Point
  {
    const CGALPoint cgal_point(1.0, 2.0, 3.0);
    const auto      dealii_point = cgal_point_to_dealii_point<3>(cgal_point);
    deallog << "deal.II Point: " << dealii_point << std::endl;
  }
}
