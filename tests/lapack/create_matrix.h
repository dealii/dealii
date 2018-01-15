// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

// auxiliary function to create a SPD FullMatrix

#include "../tests.h"

template <typename FullMatrix>
void create_spd (FullMatrix &A)
{
  const unsigned int size = A.n();
  Assert (size == A.m(), ExcDimensionMismatch(size,A.m()));

  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = i; j < size; ++j)
      {
        const double val = random_value<double>();
        Assert (val >= 0. && val <= 1.,
                ExcInternalError());
        if (i==j)
          // since A(i,j) < 1 and
          // a symmetric diagonally dominant matrix is SPD
          A(i,j) = val + size;
        else
          {
            A(i,j) = val;
            A(j,i) = val;
          }
      }
}



template <typename FullMatrix>
void create_random (FullMatrix &A)
{
  for (unsigned int i = 0; i < A.m(); ++i)
    for (unsigned int j = 0; j < A.n(); ++j)
      A(i,j) = random_value<double>();
}
