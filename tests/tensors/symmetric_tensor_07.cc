// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// in symmetric_tensor_06 we have established that contracting with a
// symmetric tensor by hand works as with a full tensor that is stored
// in non-symmetric form. here make sure that we can abbreviate the contraction

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include "../tests.h"



template <int dim>
void
test()
{
  const double            lambda = 7, mu = 5;
  SymmetricTensor<4, dim> ts;
  Tensor<4, dim>          ta;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      {
        ta[i][j][i][j] += mu;
        ta[i][j][j][i] += mu;
        ta[i][i][j][j] += lambda;
      }
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = 0; l < dim; ++l)
          ts[i][j][k][l] = ta[i][j][k][l];

  SymmetricTensor<2, dim> as, bs;
  Tensor<2, dim>          aa, ba;

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      as[i][j] = aa[i][j] = (1. + (i + 1) * (j + 1));

  bs = ts * as;
  // contract indices 2 <-> 0, 3 <-> 1
  ba = double_contract<2, 0, 3, 1>(ta, aa);

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      {
        AssertThrow(as[i][j] == aa[i][j], ExcInternalError());
        AssertThrow(bs[i][j] == ba[i][j], ExcInternalError());

        deallog << as[i][j] << ' ' << bs[i][j] << std::endl;
      }

  // test distributivity of
  // multiplication
  AssertThrow((as * ts) * as == as * (ts * as), ExcInternalError());


  // also test that the elasticity
  // tensor is positive definite
  deallog << as * ts * as << std::endl;
  Assert(as * ts * as > 0, ExcInternalError());
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
