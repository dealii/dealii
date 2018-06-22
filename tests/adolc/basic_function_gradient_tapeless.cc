// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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

// Test of basic functionality:
//  - Tapeless doubles
//  - Function and gradient computations
// Adapted from https://github.com/Homebrew/homebrew-science/blob/master/adol-c.rb
// See Adol-C manual v2.6 p65 figure 8

#include <adolc/adtl.h>
#include <adolc/drivers/drivers.h>
#include <math.h>

#include "../tests.h"

int
main(void)
{
  initlog();

  const unsigned int n = 10;
  adtl::setNumDir(n);

  adtl::adouble *x = new adtl::adouble[n];
  for (unsigned int i = 0; i < n; i++)
    {
      x[i] = (i + 1.0) / (2.0 + i);
      x[i].setADValue(i, 1);
    }

  adtl::adouble y = 1.0;
  for (unsigned int i = 0; i < n; i++)
    y *= x[i];

  // --- Function ---

  const double error_func = y.getValue() - 1.0 / (1.0 + n);

  deallog << "Error (function): " << error_func << std::endl;

  // --- Gradient ---

  double err_grad = 0;
  for (unsigned int i = 0; i < n; i++)
    err_grad += std::abs(y.getADValue(i) - y.getValue() / x[i].getValue());

  deallog << "Error (gradient): " << err_grad << std::endl;

  // -- Cleanup ---

  delete[] x;

  return 0;
}
