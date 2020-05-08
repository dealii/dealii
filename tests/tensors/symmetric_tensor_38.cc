// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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


// Check operator* with a VectorizedArray

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(5);

  const int                                        dim = 3;
  SymmetricTensor<2, dim, VectorizedArray<float>>  s1;
  SymmetricTensor<2, dim, VectorizedArray<double>> s2;

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      {
        s1[i][j] = 1 + j + i * dim;
        s2[i][j] = 2 + (j + i) * dim;
      }

  float                                            factor_float  = 2.f;
  double                                           factor_double = 3.;
  SymmetricTensor<2, dim, VectorizedArray<float>>  r1 = factor_float * s1;
  SymmetricTensor<2, dim, VectorizedArray<double>> r2 = factor_double * s2;

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      {
        deallog << i << "\t" << j << std::endl;
        deallog << r1[i][j][0] << std::endl;
        AssertThrow(std::abs(r1[i][j][0] - factor_float * s1[i][j][0]) < 1.e-10,
                    ExcInternalError());
        for (unsigned int k = 1; k < VectorizedArray<float>::size(); ++k)
          {
            AssertThrow(std::abs(r1[i][j][k] - r1[i][j][0]) < 1.e-10,
                        ExcInternalError());
          }
        deallog << r2[i][j][0] << std::endl;
        Assert(std::abs(r2[i][j][0] - factor_double * s2[i][j][0]) < 1.e-10,
               ExcInternalError());
        for (unsigned int k = 1; k < VectorizedArray<double>::size(); ++k)
          {
            AssertThrow(std::abs(r2[i][j][k] - r2[i][j][0]) < 1.e-10,
                        ExcInternalError());
          }
        deallog << std::endl;
      }

  deallog << "OK" << std::endl;
}
