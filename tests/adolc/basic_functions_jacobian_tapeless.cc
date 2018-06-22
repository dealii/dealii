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
//  - Multiple dependent functions and Jacobian computations
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

  const unsigned int m = 5;  // Dependents
  const unsigned int n = 10; // Independents
  adtl::setNumDir(n);

  adtl::adouble *x = new adtl::adouble[n];
  for (unsigned int i = 0; i < n; i++)
    {
      x[i] = (i + 1.0) / (2.0 + i);
      x[i].setADValue(i, 1);
    }

  adtl::adouble *y = new adtl::adouble[m];
  for (unsigned int j = 0; j < m; ++j)
    y[j] = 1.0;

  for (unsigned int i = 0; i < n; i++)
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
  for (unsigned int j = 0; j < m; j++)
    {
      for (unsigned int i = 0; i < n; i++)
        deallog << y[j].getADValue(i) << (i < n - 1 ? "," : "");

      deallog << std::endl;
    }

  // -- Cleanup ---

  delete[] x;
  delete[] y;

  return 0;
}
