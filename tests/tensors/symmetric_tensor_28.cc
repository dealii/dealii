// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test construction, indexing, and conversion of SymmetricTensor<2,dim>
// from/to Tensor<2,dim> with dim > 3

#include <deal.II/base/symmetric_tensor.h>

#include "../tests.h"


template <int dim>
void
test()
{
  deallog << "dim=" << dim << std::endl;

  SymmetricTensor<2, dim> s;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      s[i][j] = (i + 1) * (j + 1);

  Tensor<2, dim>          t = s;
  SymmetricTensor<2, dim> u(t);

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      deallog << (i + 1) * (j + 1) << ' ' << (int)s[i][j] << ' ' << (int)t[i][j]
              << ' ' << (int)u[i][j] << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test<3>();
  test<5>();

  deallog << "OK" << std::endl;
}
