// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test of basic functionality:
//  - Tapeless doubles
//  - Function and gradient computations
// Adapted from
// https://github.com/Homebrew/homebrew-science/blob/master/adol-c.rb See ADOL-C
// manual v2.6 p65 figure 8

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
  for (unsigned int i = 0; i < n; ++i)
    {
      x[i] = (i + 1.0) / (2.0 + i);
      x[i].setADValue(i, 1);
    }

  adtl::adouble y = 1.0;
  for (unsigned int i = 0; i < n; ++i)
    y *= x[i];

  // --- Function ---

  const double error_func = y.getValue() - 1.0 / (1.0 + n);

  deallog << "Error (function): " << error_func << std::endl;

  // --- Gradient ---

  double err_grad = 0;
  for (unsigned int i = 0; i < n; ++i)
    err_grad += std::abs(y.getADValue(i) - y.getValue() / x[i].getValue());

  deallog << "Error (gradient): " << err_grad << std::endl;

  // -- Cleanup ---

  delete[] x;

  return 0;
}
