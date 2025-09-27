// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Tests compute_eigenvalues() and eigenvalue() of TridiagonalMatrix

#include <deal.II/lac/tridiagonal_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>

#include "../tests.h"


// Build the one dimensional discrete Laplacian for n intervals
template <typename number>
void
test_laplacian(unsigned int n)
{
  TridiagonalMatrix<number> M(n - 1, true);
  for (unsigned int i = 0; i < n - 2; ++i)
    {
      M(i, i)     = 2.;
      M(i, i + 1) = -1.;
    }
  M(n - 2, n - 2) = 2.;

  M.compute_eigenvalues();
  for (unsigned int i = 0; i < 5; ++i)
    deallog << '\t' << M.eigenvalue(i) * n * n;
  deallog << "\t cond " << M.eigenvalue(n - 2) / M.eigenvalue(0) << std::endl;
}


int
main()
{
  initlog();

  test_laplacian<double>(10);
  test_laplacian<double>(20);
  test_laplacian<double>(40);
  test_laplacian<double>(80);

  test_laplacian<float>(10);
  test_laplacian<float>(20);
  test_laplacian<float>(40);
  test_laplacian<float>(80);
}
