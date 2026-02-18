// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// check serialization for Tensor<1,dim>

#include <deal.II/base/point.h>

#include "serialization.h"


void
test()
{
  const unsigned int dim = 3;

  Point<dim> p1(1., 2., 3.);

  Point<dim> p2(4., 5., 6.);

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
