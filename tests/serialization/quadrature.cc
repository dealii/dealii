// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check serialization for Quadrature

#include <deal.II/base/quadrature.h>

#include <boost/serialization/vector.hpp>

#include "serialization.h"

void
test()
{
  const unsigned int dim = 2;

  std::vector<Point<dim>> points1;
  points1.push_back(Point<dim>(0, 0));
  points1.push_back(Point<dim>(0, 1));
  points1.push_back(Point<dim>(1, 1));
  points1.push_back(Point<dim>(1, 0));

  double              w1[4] = {0.20, 0.30, 0.15, 0.35};
  std::vector<double> weights1(w1, &w1[4]);

  Quadrature<dim> q1(points1, weights1);

  std::vector<Point<dim>> points2;
  points2.push_back(Point<dim>(1, 1));
  points2.push_back(Point<dim>(1, 0));
  points2.push_back(Point<dim>(0, 0));
  points2.push_back(Point<dim>(0, 1));

  double              w2[4] = {0.1, 0.2, 0.3, 0.4};
  std::vector<double> weights2(w2, &w2[4]);

  Quadrature<dim> q2(points2, weights2);

  Quadrature<dim> q3;

  verify(q1, q2);

  verify(q1, q3);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
