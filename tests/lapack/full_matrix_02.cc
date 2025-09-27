// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// LAPACKFullMatrix::compute_lu_factorization
// LAPACKFullMatrix::solve

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>

#include "../tests.h"


// Fill a matrix with the values of the Hilbert matrix
template <typename number>
void
hilbert(LAPACKFullMatrix<number> &M, const bool nonsymmetric)
{
  const unsigned int n = M.n_rows();
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
      M(i, j) = 1. / (i + j + 1.) * (nonsymmetric ? ((i < j) ? -1. : 1.) : 1.);
}


// Multiply some vectors with the matrix and its transpose, then
// compute and apply LU factorization and see if the results are equal
// to the original vector.

void
test(const unsigned int size, const bool nonsymmetric)
{
  LAPACKFullMatrix<double> M(size, size);
  hilbert(M, nonsymmetric);

  Vector<double> u(size);
  Vector<double> v(size);
  Vector<double> x(size);
  Vector<double> y(size);

  for (unsigned int i = 0; i < size; ++i)
    {
      u(i) = i + 2.;
      x(i) = i + 2.;
    }
  M.vmult(v, u);
  M.Tvmult(y, x);
  M.compute_lu_factorization();
  M.solve(v, false);
  M.solve(y, true);

  v -= u;
  y -= x;

  deallog << v.l2_norm() << std::endl;
  deallog << y.l2_norm() << std::endl;
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  test(4, false);
  test(4, true);
  test(5, false);
  test(5, true);
}
