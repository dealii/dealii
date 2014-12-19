// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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


// check serialization for Quadrature

#include "serialization.h"
#include <deal.II/base/quadrature.h>
#include <boost/serialization/vector.hpp>

void test ()
{
  const unsigned int dim = 2;

  std::vector<Point <dim> > points1;
  points1.push_back(Point<dim>(0,0));
  points1.push_back(Point<dim>(0,1));
  points1.push_back(Point<dim>(1,1));
  points1.push_back(Point<dim>(1,0));

  double w1[4] = {0.20, 0.30, 0.15, 0.35};
  std::vector<double> weights1(w1, &w1[4]);

  Quadrature<dim> q1(points1, weights1);

  std::vector<Point <dim> > points2;
  points2.push_back(Point<dim>(1,1));
  points2.push_back(Point<dim>(1,0));
  points2.push_back(Point<dim>(0,0));
  points2.push_back(Point<dim>(0,1));

  double w2[4] = {0.1, 0.2, 0.3, 0.4};
  std::vector<double> weights2(w2, &w2[4]);

  Quadrature<dim> q2(points2, weights2);

  Quadrature<dim> q3;

  verify (q1, q2);

  verify (q1, q3);
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
