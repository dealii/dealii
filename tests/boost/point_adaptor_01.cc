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

// Check that boost adaptors for points work

#include <deal.II/boost_adaptors/point.h>

#include <boost/geometry/algorithms/equals.hpp>
#include <boost/geometry/algorithms/make.hpp>

#include "../tests.h"

namespace bg = boost::geometry;

using bgPoint = bg::model::point<double, 3, bg::cs::cartesian>;

int
main(int argc, char **argv)
{
  initlog();

  Point<3> p(1, 2, 3);
  if (bg::equals(p, bg::make<Point<3>>(1, 2, 3)) == false)
    deallog << "NOT OK" << std::endl;
  else
    deallog << "OK" << std::endl;

  if (bg::equals(p, bg::make<bgPoint>(1., 2., 3.)) == false)
    deallog << "NOT OK" << std::endl;
  else
    deallog << "OK" << std::endl;
}
