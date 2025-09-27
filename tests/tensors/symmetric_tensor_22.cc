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
    const unsigned int      dim = 2;
    SymmetricTensor<2, dim> t;
    t[0][0] = 1;
    t[1][1] = 2;
    t[0][1] = 3;

    deallog << t << std::endl;
  }

  {
    const unsigned int      dim = 3;
    SymmetricTensor<4, dim> t;
    t[0][0][0][0] = t[1][0][1][0] = t[1][1][1][1] = t[2][2][2][2] =
      t[2][0][2][0]                               = 3;

    deallog << t << std::endl;
  }

  deallog << "OK" << std::endl;
}
