// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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
  for (unsigned int i = 0; i < N; i++)
    {
      x[i] = 0.1 * i;
      y[i] = k_in * x[i] + b_in;
    }


  std::pair<double, double> fit = FESeries::linear_regression(x, y);

  deallog << "exact:      " << k_in << " " << b_in << std::endl;
  deallog << "calculated: " << fit.first << " " << fit.second << std::endl;
}



int
main()
{
  initlog();

  test();
}
