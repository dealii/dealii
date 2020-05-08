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


// Tests LAPACKFullMatrix::solve in two different variants

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>

#include "../tests.h"



void
test()
{
  const unsigned int       n = 11;
  LAPACKFullMatrix<double> A(n, n);
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
      A(i, j) = random_value<double>();
  A.compute_lu_factorization();
  LAPACKFullMatrix<double> rhs_orig(n, 3);
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < n; ++j)
      rhs_orig(j, i) = random_value<double>(-0.9923, 1.3277);

  for (unsigned int transpose = 0; transpose < 2; ++transpose)
    {
      LAPACKFullMatrix<double> rhs(rhs_orig);
      A.solve(rhs, transpose);
      for (unsigned int i = 0; i < 3; ++i)
        {
          Vector<double> check(n);
          for (unsigned int j = 0; j < n; ++j)
            check(j) = rhs_orig(j, i);
          A.solve(check, transpose);
          for (unsigned int j = 0; j < n; ++j)
            Assert(std::abs(check(j) - rhs(j, i)) < 1e-13, ExcInternalError());
        }
    }

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  test();
}
