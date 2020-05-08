// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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


// Tests LAPACKFullMatrix::operator*= and operator/=

#include <deal.II/lac/lapack_full_matrix.h>

#include <iostream>

#include "../tests.h"



void
test(const bool is_singular)
{
  const unsigned int       n = 10;
  LAPACKFullMatrix<double> A(n, n);
  A = 0;

  // generate a singular matrix if is_singular == true
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
      if (i == j && (i != n - 1 || !is_singular))
        A(i, j) = 1.0;

  try
    {
      A.compute_lu_factorization();
    }
  catch (std::exception &exc)
    {
      deallog << "matrix is singular" << std::endl;
      // Some implementations of LAPACK do not detect that the vector we pass
      // down to the compute_factorization call is compatible with a singular
      // matrix and divide by zero that triggers a floating point exception,
      // so we add a small number to the diagonal of the matrix.
      A(n - 1, n - 1) += 1e-50;
    }

  Vector<double> v(n);
  v        = 1.0;
  v(n - 1) = 0.0;
  // LU factorization can only be applied if state == lu!
  A.solve(v, false);

  deallog << "apply lu factorization succeeded with norm " << v.l2_norm()
          << std::endl;
}

int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  test(false);
  test(true);
}
