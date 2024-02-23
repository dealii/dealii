// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for Point::operator()

#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

namespace bg = boost::geometry;

void
check3d()
{
  bg::model::point<double, 3, boost::geometry::cs::cartesian> point(-1.0,
                                                                    2.0,
                                                                    0.15);

  Point<3> p(point);

  for (unsigned int i = 0; i < 3; ++i)
    deallog << p[i] << ' ';
  deallog << std::endl;
}

void
check2d()
{
  bg::model::point<double, 2, bg::cs::cartesian> point(2.0, 3.0);

  Point<2> p(point);

  for (unsigned int i = 0; i < 2; ++i)
    deallog << p[i] << ' ';
  deallog << std::endl;
}

void
check1d()
{
  bg::model::point<double, 1, bg::cs::cartesian> point(12.0);

  Point<1> p(point);

  deallog << p[0] << ' ';
  deallog << std::endl;
}

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check1d();
  check2d();
  check3d();
}
