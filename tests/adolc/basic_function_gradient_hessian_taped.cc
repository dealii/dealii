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
//  - Taped doubles
//  - Gradient and hessian computations
// Adapted from https://github.com/Homebrew/homebrew-science/blob/master/adol-c.rb

#include <adolc/adouble.h>
#include <adolc/drivers/drivers.h>
#include <adolc/taping.h>
#include <math.h>

#include "../tests.h"

int
main(void)
{
  initlog();

  const unsigned int n = 10;
  std::size_t        tape_stats[STAT_SIZE];

  double * xp = new double[n];
  double   yp = 0.0;
  adouble *x  = new adouble[n];
  adouble  y  = 1.0;

  for (unsigned int i = 0; i < n; i++)
    xp[i] = (i + 1.0) / (2.0 + i);

  trace_on(1);
  for (unsigned int i = 0; i < n; i++)
    {
      x[i] <<= xp[i];
      y *= x[i];
    }
  y >>= yp;
  delete[] x;
  trace_off();
  tapestats(1, tape_stats);

  // --- Function ---
  double *f = new double;
  function(1, 1, n, xp, f);

  const double error_func_1 = yp - 1.0 / (1.0 + n);
  const double error_func_2 = *f - 1.0 / (1.0 + n);

  deallog << "Error (function 1): " << error_func_1 << std::endl;
  deallog << "Error (function 2): " << error_func_2 << std::endl;

  // --- Gradient ---

  double *g = new double[n];
  gradient(1, n, xp, g);

  double err_grad = 0;
  for (unsigned int i = 0; i < n; i++)
    err_grad += std::abs(g[i] - yp / xp[i]);

  deallog << "Error (gradient): " << err_grad << std::endl;

  // --- Hessian ---

  double **H = new double *[n];
  for (unsigned int i = 0; i < n; ++i)
    H[i] = new double[i + 1];

  hessian(1, n, xp, H);

  double error_hess = 0;
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < n; j++)
      if (i > j)
        error_hess += std::abs(H[i][j] - g[i] / xp[j]);

  deallog << "Error (hessian): " << error_hess << std::endl;

  // -- Cleanup ---

  delete f;
  f = nullptr;

  delete[] g;
  g = nullptr;

  for (unsigned int i = 0; i < n; i++)
    delete[] H[i];
  delete[] H;
  H = nullptr;

  delete[] xp;

  return 0;
}
