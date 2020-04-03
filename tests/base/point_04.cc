// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2019 by the deal.II authors
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
    deallog << p(i) << ' ';
  deallog << std::endl;
}

void
check2d()
{
  bg::model::point<double, 2, bg::cs::cartesian> point(2.0, 3.0);

  Point<2> p(point);

  for (unsigned int i = 0; i < 2; ++i)
    deallog << p(i) << ' ';
  deallog << std::endl;
}

void
check1d()
{
  bg::model::point<double, 1, bg::cs::cartesian> point(12.0);

  Point<1> p(point);

  deallog << p(0) << ' ';
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
