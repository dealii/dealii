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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Test of basic functionality:
//  - Taped doubles
//  - Multiple dependent functions and Jacobian computations
//  - Update to using C++ data structures where possible
//  - Test that the tracing of active dependent/independent variables still
//    works with these structures
// Adapted from https://github.com/Homebrew/homebrew-science/blob/master/adol-c.rb

#include "../tests.h"
#include <adolc/adouble.h>
#include <adolc/drivers/drivers.h>
#include <adolc/taping.h>

#include <math.h>

void
test_reset_vector_values (const bool reset_values, const int tape_index)
{
  const unsigned int m = 5;  // Dependents
  const unsigned int n = 10; // Independents

  std::vector<double> xp (n, 0.0);
  std::vector<double> yp (m, 0.0);
  std::vector<adouble> x (n, 1.0); // Dep. variable values initially set here
  std::vector<adouble> y (m, 1.0);

  if (reset_values == false)
    for (unsigned int i = 0; i < n; i++)
      xp[i] = (i + 1.0) / (2.0 + i);

  trace_on(tape_index);
  for (unsigned int i = 0; i < n; i++)
    {
      x[i] <<= xp[i];
      for (unsigned int j = 0; j < m; ++j)
        y[j] *= (j+1)*x[i];
    }
  for (unsigned int j = 0; j < m; ++j)
    y[j] >>= yp[j];

  trace_off();
  // tapestats(1, tape_stats);

  // --- Change values ---
  if (reset_values == true)
    for (unsigned int i = 0; i < n; i++)
      xp[i] = (i + 1.0) / (2.0 + i);

  // --- Functions ---
  double *f = new double[m];
  function(tape_index, m, n, xp.data(), f);

  deallog << "Evaluation points:" << std::endl;
  for (unsigned int i = 0; i < n; ++i)
    deallog
        << "  x[" << i << "]: " << xp[i]
        << std::endl;

  deallog << "Function values:" << std::endl;
  for (unsigned int j = 0; j < m; ++j)
    deallog
        << "  f[" << j << "]: " << f[j]
        << "  y[" << j << "]: " << yp[j]
        << std::endl;

  // --- Jacobian ---

  double **J = new double*[m];
  for (unsigned int j = 0; j < m; ++j)
    J[j] = new double[n];

  jacobian(tape_index, m, n, xp.data(), J);

  deallog << "Function jacobian J:" << std::endl;
  for (unsigned int j = 0; j < m; j++)
    {
      for (unsigned int i = 0; i < n; i++)
        deallog << J[j][i] << (i<n-1?",":"");

      deallog << std::endl;
    }

  // -- Cleanup ---

  delete[] f;
  f = nullptr;

  for (unsigned int j = 0; j < m; j++)
    delete[] J[j];
  delete[] J;
  J = nullptr;
}

int
main(void)
{
  initlog();

  deallog << "Setting dependent variables a priori" << std::endl;
  test_reset_vector_values(false, 1); // This works
  deallog << std::endl;
  deallog << "Setting dependent variables after function definition" << std::endl;
  test_reset_vector_values(true, 2); // This doesn't
}
