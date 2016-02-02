// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2016 by the deal.II authors
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


//check method mTmult of FullMatrix on larger size than full_matrix_02 where
//we interface to the external BLAS

#include "../tests.h"

#include <deal.II/lac/full_matrix.h>


template <typename Number>
void test()
{
  FullMatrix<Number> A(2, 76), B(3, 76), C(2, 3), D(2, 3);
  for (unsigned int i=0; i<A.m(); ++i)
    for (unsigned int j=0; j<A.n(); ++j)
      A(i,j) = (double)Testing::rand()/RAND_MAX;
  for (unsigned int i=0; i<B.m(); ++i)
    for (unsigned int j=0; j<B.n(); ++j)
      B(i,j) = (double)Testing::rand()/RAND_MAX;

  A.mTmult(C, B); // C = A * B^T

  for (unsigned int i=0; i<A.m(); ++i)
    for (unsigned int j=0; j<B.m(); ++j)
      for (unsigned int k=0; k<A.n(); ++k)
        D(i,j) += A(i,k) * B(j,k);
  C.add(-1., D);

  deallog << "Difference: " << (double)C.l1_norm() << std::endl;
}

int main()
{
  initlog();
  deallog.threshold_double(1e-5);

  test<double>();
  test<float>();
  test<long double>();
}
