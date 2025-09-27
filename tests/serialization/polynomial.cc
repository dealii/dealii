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


// check serialization for Polynomial

#include <deal.II/base/polynomial.h>

#include <boost/serialization/vector.hpp>

#include "serialization.h"


void
test()
{
  double              c1[3] = {1., 2., 3.};
  std::vector<double> coefficients1(c1, &c1[3]);

  Polynomials::Polynomial<double> p1(coefficients1);

  double              c2[3] = {4., 5., 6.};
  std::vector<double> coefficients2(c2, &c2[3]);

  Polynomials::Polynomial<double> p2(coefficients2);

  Polynomials::Polynomial<double> p3;

  verify(p1, p2);

  verify(p1, p3);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
