// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2018 by the deal.II authors
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


template <int dim>
void
check()
{
  static_assert(
    Point<dim>((Point<dim>::unit_vector(0) + Point<dim>::unit_vector(0)) -
               2. * Point<dim>::unit_vector(0))(0) == 0.);
  Point<dim> p;
  for (unsigned int i = 0; i < dim; ++i)
    p[i] = i;

  for (unsigned int i = 0; i < dim; ++i)
    deallog << p(i) << ' ';
  deallog << std::endl;
}

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  static_assert((Point<1>() / 1.).norm_square() == 0.);
  static_assert(Point<2>{0., 1.}.distance_square(Point<2>{}) == 1.);
  static_assert(Point<3>(0., 0., 1.).square() == 1.);

  check<1>();
  check<2>();
  check<3>();
}
