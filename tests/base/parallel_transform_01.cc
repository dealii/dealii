// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test parallel::transform

#include <deal.II/base/parallel.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



int
main()
{
  initlog();

  const unsigned int N = 10000;
  Vector<double>     x(N), y(N);

  for (unsigned int i = 0; i < N; ++i)
    x(i) = i;

  // set y=2*x
  parallel::transform(
    x.begin(), x.end(), y.begin(), [](double i) { return 2. * i; }, 10);

  // compute y=0 from the previous result
  y -= x;
  y -= x;

  AssertThrow(y.l2_norm() == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}
