// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that boost adaptors for bounding boxes work

#include <deal.II/boost_adaptors/bounding_box.h>

#include <boost/geometry/algorithms/equals.hpp>

#include "../tests.h"

namespace bg = boost::geometry;

using bgPoint = bg::model::point<double, 3, bg::cs::cartesian>;
using bgBox   = bg::model::box<bgPoint>;

int
main(int argc, char **argv)
{
  initlog();

  Point<3> p1;
  Point<3> p2(1, 2, 3);

  BoundingBox<3> box({p1, p2});

  bgBox boost_box(bgPoint(0, 0, 0), bgPoint(1, 2, 3));

  if (bg::equals(boost_box, box) == false)
    deallog << "NOT OK" << std::endl;
  else
    deallog << "OK" << std::endl;
}
