// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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
