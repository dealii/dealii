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


// check serialization for Polynomial

#include "serialization.h"
#include <boost/serialization/vector.hpp>
#include <deal.II/base/polynomial.h>


void test ()
{
  double c1[3] = {1., 2., 3.};
  std::vector<double> coefficients1(c1, &c1[3]);

  Polynomials::Polynomial<double> p1(coefficients1);

  double c2[3] = {4., 5., 6.};
  std::vector<double> coefficients2(c2, &c2[3]);

  Polynomials::Polynomial<double> p2(coefficients2);

  Polynomials::Polynomial<double> p3;

  verify (p1, p2);

  verify (p1, p3);
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
