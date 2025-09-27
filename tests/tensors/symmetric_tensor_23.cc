// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check operator<< for SymmetricTensor<2,dim> and SymmetricTensor<4,dim>

#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  {
    SymmetricTensor<2, 1> t;
    t[0][0] = 1;

    double x[1] = {1};
    AssertThrow((t == SymmetricTensor<2, 1>(x)), ExcInternalError());
  }

  {
    SymmetricTensor<2, 2> t;
    t[0][0] = 1;
    t[1][1] = 2;
    t[0][1] = 3;

    double x[3] = {1, 2, 3};
    AssertThrow((t == SymmetricTensor<2, 2>(x)), ExcInternalError());
  }

  {
    SymmetricTensor<2, 3> t;
    t[0][0] = 1;
    t[1][1] = 2;
    t[2][2] = 3;
    t[0][1] = 4;
    t[0][2] = 5;
    t[1][2] = 6;

    double x[6] = {1, 2, 3, 4, 5, 6};
    AssertThrow((t == SymmetricTensor<2, 3>(x)), ExcInternalError());
  }

  deallog << "OK" << std::endl;
}
