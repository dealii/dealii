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


// check serialization for QIterated

#include <deal.II/base/quadrature.h>

#include <boost/serialization/vector.hpp>

#include "serialization.h"

void
test()
{
  const unsigned int dim      = 2;
  unsigned int       n_copies = 3;

  std::vector<Point<1>> points1;
  points1.push_back(Point<1>(0.));
  points1.push_back(Point<1>(1.));
  double              w1[2] = {0.5, 0.5};
  std::vector<double> weights1(w1, &w1[2]);

  std::vector<Point<1>> points2;
  points2.push_back(Point<1>(0.25));
  points2.push_back(Point<1>(0.75));
  double              w2[2] = {0.4, 0.6};
  std::vector<double> weights2(w2, &w2[2]);

  Quadrature<1> qx(points1, weights1);
  Quadrature<1> qy(points2, weights2);

  QIterated<dim> q1(qx, n_copies);

  QIterated<dim> q2(qy, n_copies);

  verify(q1, q2);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
