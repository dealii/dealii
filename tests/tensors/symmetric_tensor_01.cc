// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test symmetric 2x2 tensors

#include <deal.II/base/symmetric_tensor.h>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  SymmetricTensor<2, 2> t;
  t[0][0] = 1;
  t[1][1] = 2;
  t[0][1] = 3;

  AssertThrow(t[0][1] == t[1][0], ExcInternalError());

  // check that if a single element is
  // accessed, its transpose element gets the
  // same value
  t[1][0] = 4;
  AssertThrow(t[0][1] == 4, ExcInternalError());

  // make sure transposition doesn't change
  // anything
  AssertThrow(t == transpose(t), ExcInternalError());

  // check norm of tensor
  AssertThrow(std::fabs(t.norm() - std::sqrt(1. * 1 + 2 * 2 + 2 * 4 * 4)) <
                1e-14,
              ExcInternalError());

  // make sure norm is induced by scalar
  // product
  AssertThrow(std::fabs(t.norm() * t.norm() - t * t) < 1e-14,
              ExcInternalError());

  deallog << "OK" << std::endl;
}
