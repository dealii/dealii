// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// test parallel::transform

#include <deal.II/base/parallel.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



int
main()
{
  initlog();

  const unsigned int N = 10000;
  Vector<double>     x(N), y(N), z(N), a(N);

  for (unsigned int i = 0; i < N; ++i)
    {
      x(i) = i;
      y(i) = 2 * i;
      z(i) = 3 * i;
    }

  // set a=x+y-z, which happens to be
  // zero
  parallel::transform(
    x.begin(),
    x.end(),
    y.begin(),
    z.begin(),
    a.begin(),
    [](double i, double j, double k) { return i + j - k; },
    10);

  AssertThrow(a.l2_norm() == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}
