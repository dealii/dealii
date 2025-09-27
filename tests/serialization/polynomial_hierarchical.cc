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


// check serialization for Hierarchical

#include <deal.II/base/polynomial.h>

#include <boost/serialization/vector.hpp>

#include "serialization.h"

void
test()
{
  unsigned int              degree1 = 3;
  Polynomials::Hierarchical p1(degree1);

  unsigned int              degree2 = 7;
  Polynomials::Hierarchical p2(degree2);

  verify(p1, p2);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
