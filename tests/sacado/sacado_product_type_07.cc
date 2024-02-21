// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test scalar_product between tensors and symmetric tensors with Sacado


#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/differentiation/ad/sacado_product_types.h>

#include "../tests.h"


int
main()
{
  using Sdouble  = Sacado::Fad::DFad<double>;
  using SSdouble = Sacado::Fad::DFad<Sdouble>;
  initlog();


  // check product with Tensor<2,dim>
  Tensor<2, 2, SSdouble>        t;
  SymmetricTensor<2, 2, double> st;
  SSdouble                      a(2, 0, 7.0);
  SSdouble                      b(2, 1, 3.0);
  a.val() = Sdouble(2, 0, 7.0);
  b.val() = Sdouble(2, 1, 3.0);

  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      {
        t[i][j]  = 2. * a + i * j * b;
        st[i][j] = (1. + (i + 1) * (j * 2));
      }

  deallog << scalar_product(t, st) << std::endl;
  deallog << scalar_product(st, t) << std::endl;
}
