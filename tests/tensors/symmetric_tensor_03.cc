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


// test symmetric 2x2x2x2 tensors

#include <deal.II/base/symmetric_tensor.h>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  SymmetricTensor<4, 2> t;
  t[0][0][0][0] = 1;
  t[1][1][1][1] = 2;
  t[0][1][0][1] = 3;

  AssertThrow(t[0][1][0][1] == t[1][0][1][0], ExcInternalError());

  // check that if a single element is
  // accessed, its transpose element gets the
  // same value
  t[1][0][0][1] = 4;
  AssertThrow(t[0][1][1][0] == 4, ExcInternalError());

  // make sure transposition doesn't change
  // anything
  AssertThrow(t == transpose(t), ExcInternalError());

  // check norm of tensor
  deallog << t.norm() << std::endl;

  // make sure norm is induced by scalar
  // product
  double norm_sqr = 0;
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      for (unsigned int k = 0; k < 2; ++k)
        for (unsigned int l = 0; l < 2; ++l)
          norm_sqr += t[i][j][k][l] * t[i][j][k][l];

  AssertThrow(std::fabs(t.norm() * t.norm() - norm_sqr) < 1e-14,
              ExcInternalError());

  deallog << "OK" << std::endl;
}
