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


// check serialization for Vector

#include <deal.II/lac/vector.h>

#include <boost/serialization/vector.hpp>

#include "serialization.h"

void
test()
{
  unsigned int n = 5;

  Vector<double> v1(n);
  Vector<double> v2(n);

  Vector<double> v3;

  for (unsigned int i = 0; i < n; ++i)
    {
      v1(i) = i * 1.;
      v2(i) = i * 1. + n * 1.;
    }


  verify(v1, v2);

  verify(v1, v3);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
