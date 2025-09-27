// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test full 2x2 tensors

#include <deal.II/base/tensor.h>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  Tensor<2, 2> t;
  t[0][0] = 1;
  t[1][1] = 2;
  t[0][1] = 4;
  t[1][0] = 4;
  // make sure transposition doesn't change
  // anything
  AssertThrow(t == transpose(t), ExcInternalError());

  // check norm of tensor
  AssertThrow(std::fabs(t.norm() - std::sqrt(1. * 1 + 2 * 2 + 2 * 4 * 4)) <
                1e-14,
              ExcInternalError());

  deallog << "OK" << std::endl;
}
