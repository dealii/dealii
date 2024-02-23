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


// check serialization for Tensor<1,dim>

#include <deal.II/base/tensor.h>

#include "serialization.h"


void
test()
{
  const unsigned int dim = 3;

  double         a1[3] = {1, 2, 3};
  Tensor<1, dim> t1(a1);

  double         a2[3] = {3, 6, 9};
  Tensor<1, dim> t2(a2);

  verify(t1, t2);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
