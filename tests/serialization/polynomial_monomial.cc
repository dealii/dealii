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


// check serialization for Monomial

#include <deal.II/base/polynomial.h>

#include <boost/serialization/vector.hpp>

#include "serialization.h"


void
test()
{
  unsigned int n1           = 3;
  double       coefficient1 = 5.;

  Polynomials::Monomial<double> m1(n1, coefficient1);

  unsigned int n2           = 3;
  double       coefficient2 = 2.;

  Polynomials::Monomial<double> m2(n2, coefficient2);

  verify(m1, m2);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
