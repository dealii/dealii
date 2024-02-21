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


// test scalar_product between tensors and symmetric tensors

#include <deal.II/base/symmetric_tensor.h>

#include "../tests.h"


template <int dim>
void
test()
{
  deallog << "dim=" << dim << std::endl;

  // initialize symmetric and non-symmetric tensors. in the former case, we
  // only need to initialize one half
  SymmetricTensor<2, dim> s;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      s[i][j] = (1. + (i + 1) * (j * 2));

  Tensor<2, dim> t;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      t[i][j] = (1. + (i + 1) * (j * 2));

  deallog << scalar_product(s, s) << std::endl;
  AssertThrow(scalar_product(s, s) == s * s, ExcInternalError());

  deallog << scalar_product(s, t) << std::endl;
  AssertThrow(scalar_product(s, t) == s * symmetrize(t), ExcInternalError());
  AssertThrow(scalar_product(s, t) == symmetrize(t) * s, ExcInternalError());

  deallog << scalar_product(t, s) << std::endl;
  AssertThrow(scalar_product(t, s) == s * symmetrize(t), ExcInternalError());
  AssertThrow(scalar_product(t, s) == symmetrize(t) * s, ExcInternalError());
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test<1>();
  test<2>();
  test<3>();
  test<4>();

  deallog << "OK" << std::endl;
}
