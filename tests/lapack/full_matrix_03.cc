// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Tests SVD of LAPACKFullMatrix by comparing to vmult of FullMatrix

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>

#include "../tests.h"

/*
 * Eigenvalues and -vectors of this system are
 * lambda = 1     v = (1, 1, 1, 1)
 * lambda = 5     v = (1,-1, 0, 0)
 * lambda = 5     v = (0, 1,-1, 0)
 * lambda = 5     v = (0, 0, 1,-1)
 */
const double symm[4][4] = {{4., -1., -1., -1.},
                           {-1., 4., -1., -1.},
                           {-1., -1., 4., -1.},
                           {-1., -1., -1., 4.}};

const double rect[3][4] = {{4., 3., 2., 1.},
                           {5., 8., 1., -2.},
                           {11., 13., -4., -5}};


void
test_rect(unsigned int m, unsigned int n, const double *values)
{
  std::ostringstream prefix;
  prefix << m << 'x' << n;
  deallog.push(prefix.str());

  FullMatrix<double>       A(m, n, values);
  LAPACKFullMatrix<double> LA(m, n);
  LA = A;
  LA.compute_svd();

  deallog << "Singular values";
  for (unsigned int i = 0; i < LA.n_rows(); ++i)
    deallog << ' ' << LA.singular_value(i);
  deallog << std::endl;

  Vector<double> u(n);
  Vector<double> v1(m);
  Vector<double> v2(m);

  for (unsigned int i = 0; i < u.size(); ++i)
    u(i) = i * i;

  // Test rectangular vmult. All
  // results compare with same
  // operation for FullMatrix.

  A.vmult(v1, u);
  LA.vmult(v2, u);
  v1 -= v2;
  if (v1.l2_norm() < 1.e-12)
    deallog << "vmult ok" << std::endl;
  else
    deallog << "vmult error " << v1.l2_norm() << std::endl;
  v1 = v2;

  A.vmult_add(v1, u);
  LA.vmult_add(v2, u);
  v1 -= v2;
  if (v1.l2_norm() < 1.e-12)
    deallog << "vmult_add ok" << std::endl;
  else
    deallog << "vmult_add error " << v1.l2_norm() << std::endl;

  LA.Tvmult(u, v2);
  u *= -1;
  A.Tvmult_add(u, v2);
  if (u.l2_norm() < 1.e-11)
    deallog << "Tvmult ok" << std::endl;
  else
    deallog << "Tvmult error " << u.l2_norm() << std::endl;

  A.Tvmult(u, v2);
  u *= -1;
  LA.Tvmult_add(u, v2);
  if (u.l2_norm() < 1.e-11)
    deallog << "Tvmult_add ok" << std::endl;
  else
    deallog << "Tvmult_add error " << u.l2_norm() << std::endl;

  deallog.pop();
}

int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  test_rect(4, 4, &symm[0][0]);
  test_rect(4, 3, &rect[0][0]);
  test_rect(3, 4, &rect[0][0]);

  // Test symmetric system
  FullMatrix<double>       A(4, 4, &symm[0][0]);
  LAPACKFullMatrix<double> LA(4, 4);
  LA = A;
  LA.compute_eigenvalues();
  for (unsigned int i = 0; i < A.m(); ++i)
    {
      std::complex<double> lambda = LA.eigenvalue(i);
      deallog << "Eigenvalues " << (int)(lambda.real() + .0001) << '\t'
              << (int)(lambda.imag() + .0001) << std::endl;
    }
}
