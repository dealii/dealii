//
// Copyright (C) 2016 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Test of basic functionality:
//  - Taped doubles
//  - Multiple dependent functions and Jacobian computations
// Adapted from https://github.com/Homebrew/homebrew-science/blob/master/adol-c.rb

#include "../tests.h"
#include <adolc/adouble.h>
#include <adolc/drivers/drivers.h>
#include <adolc/taping.h>

#include <math.h>

int
main(void)
{
  initlog();

  const unsigned int m = 5;  // Dependents
  const unsigned int n = 10; // Independents
  std::size_t        tape_stats[STAT_SIZE];

  double*  xp = new double[n];
  double*  yp = new double[m];
  adouble* x  = new adouble[n];
  adouble* y  = new adouble[m];

  for(unsigned int i = 0; i < n; i++)
    xp[i] = (i + 1.0) / (2.0 + i);

  for(unsigned int j = 0; j < m; ++j)
    y[j] = 1.0;

  trace_on(1);
  for(unsigned int i = 0; i < n; i++)
    {
      x[i] <<= xp[i];
      for(unsigned int j = 0; j < m; ++j)
        y[j] *= (j + 1) * x[i];
    }
  for(unsigned int j = 0; j < m; ++j)
    y[j] >>= yp[j];

  delete[] x;
  delete[] y;

  trace_off();
  tapestats(1, tape_stats);

  // --- Functions ---
  double* f = new double[m];
  function(1, m, n, xp, f);

  deallog << "Evaluation points:" << std::endl;
  for(unsigned int i = 0; i < n; ++i)
    deallog << "  x[" << i << "]: " << xp[i] << std::endl;

  deallog << "Function values:" << std::endl;
  for(unsigned int j = 0; j < m; ++j)
    deallog << "  f[" << j << "]: " << f[j] << "  y[" << j << "]: " << yp[j]
            << std::endl;

  // --- Jacobian ---

  double** J = new double*[m];
  for(unsigned int j = 0; j < m; ++j)
    J[j] = new double[n];

  jacobian(1, m, n, xp, J);

  deallog << "Function jacobian J:" << std::endl;
  for(unsigned int j = 0; j < m; j++)
    {
      for(unsigned int i = 0; i < n; i++)
        deallog << J[j][i] << (i < n - 1 ? "," : "");

      deallog << std::endl;
    }

  // -- Cleanup ---

  delete[] xp;
  delete[] yp;

  delete[] f;
  f = nullptr;

  for(unsigned int j = 0; j < m; j++)
    delete[] J[j];
  delete[] J;
  J = nullptr;

  return 0;
}
