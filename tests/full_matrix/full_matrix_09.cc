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


// check method TmTmult of FullMatrix on larger size than full_matrix_03 where
// we interface to the external BLAS by comparison to mmult of the transpose
// matrix

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"


template <typename Number>
void
test()
{
  FullMatrix<Number> A(2, 76), B(76, 3), C(2, 3), D(3, 2), E(2, 3);
  for (unsigned int i = 0; i < A.m(); ++i)
    for (unsigned int j = 0; j < A.n(); ++j)
      A(i, j) = random_value<double>();
  for (unsigned int i = 0; i < B.m(); ++i)
    for (unsigned int j = 0; j < B.n(); ++j)
      B(i, j) = random_value<double>();

  A.mmult(C, B);   // C = A * B
  B.TmTmult(D, A); // D = B^T * A^T

  E.copy_transposed(D); // E = D^T = C

  C.add(-1., E);

  const Number tolerance = 100 * std::numeric_limits<Number>::epsilon();
  deallog << "Difference 1: "
          << filter_out_small_numbers(C.l1_norm(), tolerance) << std::endl;

  C = 0;
  for (unsigned int i = 0; i < A.m(); ++i)
    for (unsigned int j = 0; j < B.n(); ++j)
      for (unsigned int k = 0; k < A.n(); ++k)
        C(i, j) += A(i, k) * B(k, j);
  C.add(-1., E);

  deallog << "Difference 2: "
          << filter_out_small_numbers(C.l1_norm(), tolerance) << std::endl;
}

int
main()
{
  initlog();

  test<double>();
  test<float>();
  test<long double>();
}
