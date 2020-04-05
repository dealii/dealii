// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

// test that multiplying a rank-4 (identity) symmetric tensor and
// a rank-2 symmetric tensor works correctly for VectorizedArray<double>

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/vectorization.h>

#include "../tests.h"


// Trying to fill the rank-2 ensor in the main function triggers a compiler bug
// in clang-3.7.0 and clang-3.9.1 in release mode. Hence, use a separate
// function.
template <int dim, typename Number>
void fill_tensor(
  dealii::SymmetricTensor<2, dim, dealii::VectorizedArray<Number>> &A)
{
  Number counter = 0.0;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int v = 0; v < dealii::VectorizedArray<Number>::size(); v++)
        {
          A[i][j][v] = counter;
          counter += 1.0;
        }
}

int
main()
{
  initlog();

  const unsigned int dim = 3;

  SymmetricTensor<4, dim, VectorizedArray<double>> I;
  SymmetricTensor<2, dim, VectorizedArray<double>> A;
  SymmetricTensor<2, dim, VectorizedArray<double>> B;

  // I^sym = 0.5(d_ik*d_jl + d_il*d_jk) -> I^sym : A = A^sym
  for (unsigned int i = 0; i < dim; i++)
    for (unsigned int j = 0; j < dim; j++)
      for (unsigned int k = 0; k < dim; k++)
        for (unsigned int l = 0; l < dim; l++)
          I[i][j][k][l] = ((i == k && j == l && i == l && j == k) ?
                             make_vectorized_array(1.0) :
                             ((i == k && j == l) || (i == l && j == k) ?
                                make_vectorized_array(0.5) :
                                make_vectorized_array(0.0)));

  fill_tensor(A);

  B = I * A;
  B -= A;

  // Note that you cannot use B.norm() here even with something
  // like VectorizedArray<double> frob_norm = B.norm() -> Maybe a TODO?
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v)
        if (B[i][j][v] != 0.0)
          deallog << "Not OK" << std::endl;

  deallog << "OK" << std::endl;
}
