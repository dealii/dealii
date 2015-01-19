// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// check serialization for QAnisotropic

#include "serialization.h"
#include <deal.II/base/quadrature.h>
#include <boost/serialization/vector.hpp>

void test ()
{
  const unsigned int dim = 2;

  std::vector<Point <1> > points1;
  points1.push_back(Point<1>(0.));
  points1.push_back(Point<1>(1.));
  double w1[2] = {0.5, 0.5};
  std::vector<double> weights1(w1, &w1[2]);

  std::vector<Point <1> > points2;
  points2.push_back(Point<1>(0.25));
  points2.push_back(Point<1>(0.75));
  double w2[2] = {0.4, 0.6};
  std::vector<double> weights2(w2, &w2[2]);

  Quadrature<1> qx(points1, weights1);
  Quadrature<1> qy(points2, weights2);

  QAnisotropic<dim> q1(qx, qy);

  QAnisotropic<dim> q2(qy, qx);

  verify (q1, q2);
}


int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
