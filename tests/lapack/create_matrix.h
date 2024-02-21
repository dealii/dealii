// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// auxiliary function to create a SPD FullMatrix

#include "../tests.h"

template <typename FullMatrix>
void
create_spd(FullMatrix &A)
{
  const unsigned int size = A.n();
  Assert(size == A.m(), ExcDimensionMismatch(size, A.m()));

  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = i; j < size; ++j)
      {
        const double val = random_value<typename FullMatrix::value_type>();
        Assert(val >= 0. && val <= 1., ExcInternalError());
        if (i == j)
          // since A(i,j) < 1 and
          // a symmetric diagonally dominant matrix is SPD
          A(i, j) = val + size;
        else
          {
            A(i, j) = val;
            A(j, i) = val;
          }
      }
}

// create random invertible lower triangular matrix
template <typename FullMatrix>
void
create_random_lt(FullMatrix &A)
{
  const unsigned int size = A.n();
  Assert(size == A.m(), ExcDimensionMismatch(size, A.m()));
  A = 0.;
  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j <= i; ++j)
      {
        if (i == j)
          A(i, j) = (typename FullMatrix::value_type)(1.);
        else
          A(i, j) = random_value<typename FullMatrix::value_type>();
      }
}

template <typename FullMatrix>
void
create_random(FullMatrix &A)
{
  for (unsigned int i = 0; i < A.m(); ++i)
    for (unsigned int j = 0; j < A.n(); ++j)
      A(i, j) = random_value<typename FullMatrix::value_type>();
}



template <typename NumberType>
void
create_random(Vector<NumberType> &V)
{
  for (unsigned int i = 0; i < V.size(); ++i)
    V(i) = random_value<NumberType>();
}
