// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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



#include <deal.II/base/quadrature.h>

#include "../tests.h"

/*
 * Test that it is possible to add a point to the quadrature formula using the
 * push_back function.
 */
template <int dim>
void
test_push_back()
{
  Quadrature<dim>  quadrature;
  const Point<dim> point  = Point<dim>::unit_vector(0);
  const double     weight = .5;
  quadrature.push_back(point, weight);

  AssertThrow(1 == quadrature.size(), ExcInternalError());
  deallog << "point: " << quadrature.point(0) << std::endl;
  deallog << "weight: " << quadrature.weight(0) << std::endl;
  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  test_push_back<1>();
  test_push_back<2>();
  test_push_back<3>();
}
