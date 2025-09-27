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
//  - Multiple dependent functions and Jacobian computations
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

  const unsigned int m = 5;  // Dependents
  const unsigned int n = 10; // Independents
  adtl::setNumDir(n);

  adtl::adouble *x = new adtl::adouble[n];
  for (unsigned int i = 0; i < n; ++i)
    {
      x[i] = (i + 1.0) / (2.0 + i);
      x[i].setADValue(i, 1);
    }

  adtl::adouble *y = new adtl::adouble[m];
  for (unsigned int j = 0; j < m; ++j)
    y[j] = 1.0;

  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < m; ++j)
      y[j] *= (j + 1) * x[i];

  // --- Functions ---

  deallog << "Evaluation points:" << std::endl;
  for (unsigned int i = 0; i < n; ++i)
    deallog << "  x[" << i << "]: " << x[i].getValue() << std::endl;

  deallog << "Function values:" << std::endl;
  for (unsigned int j = 0; j < m; ++j)
    deallog << "  y[" << j << "]: " << y[j].getValue() << std::endl;

  // --- Jacobian ---

  deallog << "Function jacobian J:" << std::endl;
  for (unsigned int j = 0; j < m; ++j)
    {
      for (unsigned int i = 0; i < n; ++i)
        deallog << y[j].getADValue(i) << (i < n - 1 ? "," : "");

      deallog << std::endl;
    }

  // -- Cleanup ---

  delete[] x;
  delete[] y;

  return 0;
}
