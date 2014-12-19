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

//
// Author: Andrew McBride
//
// check on the Polynomial::operator -=

#include <deal.II/base/logstream.h>
#include <deal.II/base/polynomial.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace dealii;


int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);

  //      subtract two equal polynomials up to order p_dim
  //      the result evaluated at an arbitrary point
  //      should always be zero

  //      check polynomials up to order 32
  const unsigned int p_dim = 32;

  std::vector<double> coefficients_a;

  //      arbitrary point
  double evalutation_number = 12.123;

  for (unsigned int pp = 0; pp < p_dim; pp++)
    {
      coefficients_a.push_back(pp);

      Polynomials::Polynomial<double> test_poly_a(coefficients_a);
      Polynomials::Polynomial<double> test_poly_b(coefficients_a);

      test_poly_b -= test_poly_a;

      deallog << test_poly_b.value(evalutation_number) << std::endl;
    }
}
