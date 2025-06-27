// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

//
// Author: Andrew McBride
//
// check on the Polynomial::operator -=

#include <deal.II/base/polynomial.h>

#include <iostream>
#include <vector>

#include "../tests.h"



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  //      subtract two equal polynomials up to order p_dim
  //      the result evaluated at an arbitrary point
  //      should always be zero

  //      check polynomials up to order 32
  const unsigned int p_dim = 32;

  std::vector<double> coefficients_a;

  //      arbitrary point
  double evaluation_number = 12.123;

  for (unsigned int pp = 0; pp < p_dim; ++pp)
    {
      coefficients_a.push_back(pp);

      Polynomials::Polynomial<double> test_poly_a(coefficients_a);
      Polynomials::Polynomial<double> test_poly_b(coefficients_a);

      test_poly_b -= test_poly_a;

      deallog << test_poly_b.value(evaluation_number) << std::endl;
    }
}
