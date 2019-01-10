// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2017 by the deal.II authors
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

#include <boost/geometry.hpp>

#include "../tests.h"

namespace bg = boost::geometry;

template <int dim>
void
check()
{
  bg::model::point<double, dim, bg::cs::cartesian> bg_point;

  constexpr unsigned int y_index = (spacedim < 2) ? 0 : 1;
  constexpr unsigned int z_index = (spacedim < 3) ? 0 : 2;

  bg_point.set<0>(42.0);

  if (dim >= 2)
    bg_point.set<y_index>(42.0 + dim);

  if (dim >= 3)
    bg_point.set<z_index>(42.0 + dim + 1);

  Point<dim> p(bg_point);

  for (unsigned int i = 0; i < dim; ++i)
    deallog << p(i) << ' ';
  deallog << std::endl;
}

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check<1>();
  check<2>();
  check<3>();
}
