// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test FESeries::linear_regression()

#include <deal.II/fe/fe_series.h>

#include <iostream>

#include "../tests.h"



void
test()
{
  const double        k_in = numbers::PI;
  const double        b_in = std::sqrt(2.);
  const unsigned int  N    = 10;
  std::vector<double> x(N), y(N);

  // fill the data
  for (unsigned int i = 0; i < N; ++i)
    {
      x[i] = 0.1 * i;
      y[i] = k_in * x[i] + b_in;
    }


  std::pair<double, double> fit = FESeries::linear_regression(x, y);

  deallog << "exact:      " << k_in << ' ' << b_in << std::endl;
  deallog << "calculated: " << fit.first << ' ' << fit.second << std::endl;
}



int
main()
{
  initlog();

  test();
}
